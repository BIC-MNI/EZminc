#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use POSIX qw(floor);
use strict;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $debug=0;
my @image_range;
my @ovl_range;

my @coronal;
my @sagittal;
my @axial;

my $overlay;
my $separate;
my $out_coronal;
my $out_axial;
my $out_sagittal;
my $mask;
my $scale;
my $max;
my $over;
my ($red,$green,$blue);
my $vertical;
my ($trim_x,$trim_y,$trim_z)=(0,0,0);
my $adjust;

GetOptions (    
          "verbose"          => \$verbose,
          "clobber"          => \$clobber,
          'image-range=f{2}' => \@image_range,
          'ovl-range=f{2}'   => \@ovl_range,
          "overlay=s"        => \$overlay,
          "mask=s"           => \$mask,
          "scale=n"          => \$scale,
          'separate'         => \$separate,
          'coronal=n{,}'       => \@coronal,
          'sagittal=n{,}'      => \@sagittal,
          'axial=n{,}'         => \@axial,
          'max'              => \$max,
          'over'             => \$over,
          'red=s'            => \$red,
          'green=s'          => \$green,
          'blue=s'           => \$blue,
          'vertical'         => \$vertical,
          'trim-x=n'         => \$trim_x,
          'trim-y=n'         => \$trim_y,
          'trim-z=n'         => \$trim_z,
          'adjust'           => \$adjust
          );
          
die <<END  if $#ARGV<1;
Usage: $me <scan_in> <scan_out> 
[ --verbose
  --clobber
  --image-range a b
  --ovl-range a b 
  --overlay <overlay>
  --separate
  --mask <mask for overlay>
  --scale <n>
  --coronal <n>
  --axial <n>
  --sagittal <n>
  --max - use max to combine images
  --red <mnc> - use overlay in red
  --green <mnc> - overlay in green
  --blue <mnc> -overlay in blue
  --vertical - montage vertically 
  --trim-x <amount> trim X direction
  --trim-y <amount> trim Y direction
  --trim-z <amount> trim Z direction
  --over put overlay on top of image (works for --red --green --blue for now)
  --adjust coordinate when using trim
]
END

my ($in,$out)=@ARGV;

print "$in $out\n";

@coronal=(116) if $#coronal==-1;
@sagittal=(93) if $#sagittal==-1;
@axial=(76)    if $#axial==-1;

if($separate) {
   my $i;  
  if($#coronal==0)
  {
    $out_coronal=$out."_coronal.png";
    check_file($out_coronal) unless $clobber;
  } else {

    for($i=0;$i<=$#coronal;$i+=1)
    {
      check_file("${out}_coronal_$i.png") unless $clobber;
    }
  }

  if($#axial==0)
  {
    $out_axial=$out."_axial.png";
    check_file($out_axial) unless $clobber;
  } else {

    for($i=0;$i<=$#axial;$i+=1)
    {
      check_file("${out}_axial_$i.png") unless $clobber;
    }
  }

  if($#sagittal==0)
  {
    $out_sagittal=$out."_sagittal.png";
    check_file($out_sagittal) unless $clobber;
  } else {

    for($i=0;$i<=$#sagittal;$i+=1)
    {
      check_file("${out}_sagittal_$i.png") unless $clobber;
    }
  }
} else {
  check_file($out) unless $clobber;
}
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};

if($trim_x || $trim_y || $trim_z)
{
  my ($xs,$ys,$zs)=split("\n",`mincinfo -dimlength xspace -dimlength yspace -dimlength zspace $in`);
  $xs*=1;
  $ys*=1;
  $zs*=1;
  my ($xl,$yl,$zl)=($xs-2*$trim_x, $ys-2*$trim_y, $zs-2*$trim_z);

  do_cmd('mincreshape',
          '-dimrange',"xspace=$trim_x,$xl",
          '-dimrange',"yspace=$trim_y,$yl",
          '-dimrange',"zspace=$trim_z,$zl",
          $in,"$tmpdir/in.mnc");

  $in="$tmpdir/in.mnc";
  if($adjust)
  {
    my $i;

    for($i=0;$i<=$#coronal;$i+=1)
    {
      $coronal[$i]-=$trim_y;
    }

    for($i=0;$i<=$#sagittal;$i+=1)
    {
      $sagittal[$i]-=$trim_x;
    }

    for($i=0;$i<=$#axial;$i+=1)
    {
      $axial[$i]-=$trim_z;
    }
  }
}

#produce RGB images
my @args=("minclookup",$in,"$tmpdir/gray.mnc",'-gray','-byte');
push @args,'-min',$image_range[0],'-max',$image_range[1] if $#image_range>0;
do_cmd(@args);

if($overlay)
{
  do_cmd('mincresample',$overlay,'-like',$in,"$tmpdir/overlay.mnc");
  if($mask)
  {
    do_cmd('mincresample',$mask,'-like',$in,'-nearest',"$tmpdir/mask.mnc");
    $mask="$tmpdir/mask.mnc";
    do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/overlay.mnc",$mask,"$tmpdir/masked_overlay.mnc");
    do_cmd('mv',"$tmpdir/masked_overlay.mnc","$tmpdir/overlay.mnc");
  }
  
  @args=("minclookup","$tmpdir/overlay.mnc","$tmpdir/overlay_rgb.mnc",'-byte');
  push @args,'-spectral';
  
  push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
  do_cmd(@args);
  
  #mix colors
  if($max) {
    do_cmd('mincmath','-max',"$tmpdir/gray.mnc","$tmpdir/overlay_rgb.mnc","$tmpdir/gray_.mnc");
  } else {
    do_cmd('mincaverage',"$tmpdir/gray.mnc","$tmpdir/overlay_rgb.mnc","$tmpdir/gray_.mnc");
  }
  do_cmd('mv',"$tmpdir/gray_.mnc","$tmpdir/gray.mnc");
}

if($red || $green || $blue)
{
  my @channels;
  if($red)
  {
    do_cmd('mincresample',$red,'-like',$in,"$tmpdir/red.mnc");
    if($mask)
    {
      do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/red.mnc",$mask,"$tmpdir/masked_red.mnc");
      do_cmd('minccalc','-express',"A[0]>0.5&&A[1]>$ovl_range[0]?1:0",'-byte',$mask,"$tmpdir/red.mnc","$tmpdir/ovl_mask.mnc") if $over && $#ovl_range>0&& !-e "$tmpdir/ovl_mask.mnc";
      
      do_cmd('mv',"$tmpdir/masked_red.mnc","$tmpdir/red.mnc");
    } else {
      do_cmd('minccalc','-express',"A[0]>$ovl_range[0]?1:0",'-byte',"$tmpdir/red.mnc","$tmpdir/ovl_mask.mnc") if $over && $#ovl_range>0;

    }
    @args=("minclookup","$tmpdir/red.mnc","$tmpdir/overlay_red.mnc",'-byte','-lut_string','0.0 0.0 0.0 0.0;1.0 1.0 0.0 0.0');
    push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
    do_cmd(@args);
    push @channels,"$tmpdir/overlay_red.mnc";
  }
  
  if($green)
  {
    do_cmd('mincresample',$green,'-like',$in,"$tmpdir/green.mnc");
    if($mask)
    {
      do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/green.mnc",$mask,"$tmpdir/masked_green.mnc");
      do_cmd('minccalc','-express',"A[0]>0.5&&A[1]>$ovl_range[0]?1:0",'-byte',$mask,"$tmpdir/green.mnc","$tmpdir/ovl_mask.mnc") if $over && $#ovl_range>0&& !-e "$tmpdir/ovl_mask.mnc";
      
      do_cmd('mv',"$tmpdir/masked_green.mnc","$tmpdir/green.mnc");
    } else {
      do_cmd('minccalc','-express',"A[0]>$ovl_range[0]?1:0",'-byte',"$tmpdir/green.mnc","$tmpdir/ovl_mask.mnc") if $over && $#ovl_range>0 && !-e "$tmpdir/ovl_mask.mnc";
    }
    @args=("minclookup","$tmpdir/green.mnc","$tmpdir/overlay_green.mnc",'-byte','-lut_string','0.0 0.0 0.0 0.0;1.0 0.0 1.0 0.0');
    push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
    do_cmd(@args);
    push @channels,"$tmpdir/overlay_green.mnc";
  }
  
  if($blue)
  {
    do_cmd('mincresample',$blue,'-like',$in,"$tmpdir/blue.mnc");
    if($mask)
    {
      do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/blue.mnc",$mask,"$tmpdir/masked_blue.mnc");
      do_cmd('minccalc','-express',"A[0]>0.5&&A[1]>$ovl_range[0]?1:0",'-byte',$mask,"$tmpdir/blue.mnc","$tmpdir/ovl_mask.mnc") if $over && $#ovl_range>0&& !-e "$tmpdir/ovl_mask.mnc";
      
      do_cmd('mv',"$tmpdir/masked_blue.mnc","$tmpdir/blue.mnc");
    } else {
      do_cmd('minccalc','-express',"A[0]>$ovl_range[0]?1:0",'-byte',"$tmpdir/blue.mnc","$tmpdir/ovl_mask.mnc") if $over && $#ovl_range>0 && !-e "$tmpdir/ovl_mask.mnc";
    }
    @args=("minclookup","$tmpdir/blue.mnc","$tmpdir/overlay_blue.mnc",'-byte','-lut_string','0.0 0.0 0.0 0.0;1.0 0.0 0.0 1.0');
    push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
    do_cmd(@args);
    push @channels,"$tmpdir/overlay_blue.mnc";
  }
   
  
  #mix colors
  if($over && -e "$tmpdir/ovl_mask.mnc") {
    do_cmd('minclookup','-gray',"$tmpdir/ovl_mask.mnc","$tmpdir/ovl_mask_g.mnc",'-min',0,'-max',1);
    do_cmd('mincmath','-max',@channels,"$tmpdir/mix.mnc");
    do_cmd('minccalc','-express','A[2]>0.5?A[1]:A[0]',"$tmpdir/gray.mnc","$tmpdir/mix.mnc","$tmpdir/ovl_mask_g.mnc","$tmpdir/gray_.mnc");
    
  } elsif($max) {
    do_cmd('mincmath','-max',"$tmpdir/gray.mnc",@channels,"$tmpdir/gray_.mnc");
  } else {
    do_cmd('mincaverage',"$tmpdir/gray.mnc",@channels,"$tmpdir/gray_.mnc");
  }
  do_cmd('mv',"$tmpdir/gray_.mnc","$tmpdir/gray.mnc");
}

# produce individual slices

my $i;
my @pics;

for($i=0;$i<=$#coronal;$i+=1)
{
  @args=('mincpik',"$tmpdir/gray.mnc","$tmpdir/coronal_$i.miff",'-coronal','-slice',$coronal[$i]);
  push @args,'-scale',$scale if $scale;
  do_cmd(@args);
  push @pics,"$tmpdir/coronal_$i.miff";
}

for($i=0;$i<=$#axial;$i+=1)
{
  @args=('mincpik',"$tmpdir/gray.mnc","$tmpdir/axial_$i.miff",'-axial','-slice',$axial[$i]);
  push @args,'-scale',$scale if $scale;
  do_cmd(@args);
  push @pics,"$tmpdir/axial_$i.miff";
}

for($i=0;$i<=$#sagittal;$i+=1)
{
  @args=('mincpik',"$tmpdir/gray.mnc","$tmpdir/sagittal_$i.miff",'-sagittal','-slice',$sagittal[$i]);
  push @args,'-scale',$scale if $scale;
  do_cmd(@args);
  push @pics,"$tmpdir/sagittal_$i.miff";
}

if($separate)
{

  if($#coronal==0)
  {
    do_cmd("convert","$tmpdir/coronal_0.miff",$out_coronal);
  } else {
    for($i=0;$i<=$#coronal;$i+=1)
    {  
      do_cmd("convert","$tmpdir/coronal_$i.miff","${out}_coronal_$i.png");
    }
  }
  
  if($#axial==0)
  {
    do_cmd("convert","$tmpdir/axial_0.miff",$out_axial);
  } else {
    for($i=0;$i<=$#axial;$i+=1)
    {  
      do_cmd("convert","$tmpdir/axial_$i.miff","${out}_axial_$i.png");
    }
  }

  if($#sagittal==0)
  {
    do_cmd("convert","$tmpdir/sagittal_0.miff",$out_sagittal);
  } else {
    for($i=0;$i<=$#sagittal;$i+=1)
    {  
      do_cmd("convert","$tmpdir/sagittal_$i.miff","${out}_sagittal_$i.png");
    }
  }

} else {
  my $geo='+0+0';

  if($vertical)
  {
    equalize_width(@pics);
  } else {
    equalize_height(@pics);
  }

  my $lines=floor(($#pics+2)/3);
  #chomp($geo);
  do_cmd('montage','-tile',$vertical?"${lines}x3":"3x$lines",
          '-geometry',$geo,'-gravity','North',
          '-background','black',
         @pics,$out);
}

#####################################################################
sub do_cmd { 
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
      system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  if($_[0] && -e $_[0])
  {
    die("${_[0]} exists!\n");
    return 0;
  }    
  return 1;
}

sub equalize_height {
  my @in=@_;
  my @h;    
  my $i;

  for($i=0;$i<=$#in;$i+=1)
  {
    chomp($h[$i]=`identify -format '%h' $in[$i] `);
  }
  
  my @hh = sort {$b <=> $a} @h;
  my $mh=$hh[0];
  
  for($i=0;$i<=$#in;$i+=1)
  {
    my $w=($mh-$h[$i])/2;
    next if $w<=0;
    do_cmd('convert','-bordercolor','black','-border',"1x${w}",$in[$i],$in[$i]);
  }
}

sub equalize_width {
  my @in=@_;
  my @h;    
  my $i;

  for($i=0;$i<=$#in;$i+=1)
  {
    chomp($h[$i]=`identify -format '%w' $in[$i] `);
  }
  
  my @hh = sort {$b <=> $a} @h;
  my $mh=$hh[0];
  
  for($i=0;$i<=$#in;$i+=1)
  {
    my $w=($mh-$h[$i])/2;
    next if $w<=0;
    do_cmd('convert','-bordercolor','black','-border',"${w}x1",$in[$i],$in[$i]);
  }
}
