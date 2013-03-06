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
my $coronal=116;
my $sagittal=93;
my $axial=76;
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
my ($shift_x,$shift_y,$shift_z)=(0,0,0);
my $adjust;
my $cyanred=0;
my $lut;
my $discrete;
my $background='black';
my $foreground='white';
my $title;
my $pointsize;

GetOptions (    
          "verbose"          => \$verbose,
          "clobber"          => \$clobber,
          'image-range=f{2}' => \@image_range,
          'ovl-range=f{2}'   => \@ovl_range,
          "overlay=s"        => \$overlay,
          "mask=s"           => \$mask,
          "scale=n"          => \$scale,
          'separate'         => \$separate,
          'coronal=n'        => \$coronal,
          'sagittal=n'       => \$sagittal,
          'axial=n'          => \$axial,
          'max'              => \$max,
          'over'             => \$over,
          'red=s'            => \$red,
          'green=s'          => \$green,
          'blue=s'           => \$blue,
          'vertical'         => \$vertical,
          'trim-x=n'         => \$trim_x,
          'trim-y=n'         => \$trim_y,
          'trim-z=n'         => \$trim_z,
          'adjust'           => \$adjust,
          'shift-x=n'        => \$shift_x,
          'shift-y=n'        => \$shift_y,
          'shift-z=n'        => \$shift_z,
          'cyanred'          => \$cyanred,
          'lut=s'            => \$lut,
          'discrete'         => \$discrete,
          'background=s'     => \$background,
          'foreground=s'     => \$foreground,
          'title=s'          => \$title,
          'pointsize=n'      => \$pointsize,
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
  --shift-x <amount> shift X direction
  --shift-y <amount> shift Y direction
  --shift-z <amount> shift Z direction
  --over put overlay on top of image (works for --red --green --blue for now)
  --adjust coordinate when using trim+shift
  --cyanred use cyan-red for overlay
  --lut <lut table for overlay>
  --discrete <use discrete labels for overlay>
  --background <color>, default black
  --foreground <color>, default white
  --tile <title>
  --pointsize <n>
]
END

my ($in,$out)=@ARGV;

if($separate) {
  $out_coronal=$out."_coronal.png";
  $out_axial=$out."_axial.png";
  $out_sagittal=$out."_sagittal.png";

  check_file($out_coronal) unless $clobber;
  check_file($out_axial) unless $clobber;
  check_file($out_sagittal) unless $clobber;
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
  
  $trim_x+=$shift_x;
  $trim_y+=$shift_y;
  $trim_z+=$shift_z;
  
  do_cmd('mincreshape',
          '-dimrange',"xspace=$trim_x,$xl",
          '-dimrange',"yspace=$trim_y,$yl",
          '-dimrange',"zspace=$trim_z,$zl",
          $in,"$tmpdir/in.mnc");

  $in="$tmpdir/in.mnc";
  if($adjust)
  {
      $coronal-=$trim_y;
      $sagittal-=$trim_x;
      $axial-=$trim_z;
  }
}

#produce RGB images
my @args=("minclookup",$in,"$tmpdir/gray.mnc",'-gray','-byte');
push @args,'-min',$image_range[0],'-max',$image_range[1] if $#image_range>0;
do_cmd(@args);

if($overlay)
{
#   my $rr=`mincstats -q -min -max $overlay`;
#   my @range=split("\n",$rr);
  
  do_cmd('mincresample',$overlay,'-like',$in,"$tmpdir/overlay.mnc",'-float');
  if($mask)
  {
    do_cmd('mincresample',$mask,'-like',$in,'-nearest',"$tmpdir/mask.mnc");
    $mask="$tmpdir/mask.mnc";
    do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/overlay.mnc",$mask,"$tmpdir/masked_overlay.mnc");
    do_cmd('mv',"$tmpdir/masked_overlay.mnc","$tmpdir/overlay_.mnc");
  }
  
  @args=("minclookup","$tmpdir/overlay.mnc","$tmpdir/overlay_rgb.mnc",'-byte');

  if($cyanred)
  { 
    push @args,'-lut_string',"0.000 0.8 1.0 1.0;0.125 0.4 0.9 1.0;0.250 0.0 0.6 1.0;\
0.375 0.0 0.2 0.5;0.500 0.0 0.0 0.0;0.625 0.5 0.0 0.0;0.750 1.0 0.4 0.0;\
0.825 1.0 0.8 0.4;1.000 1.0 0.8 0.8";

  } elsif($lut) {
    
    push @args,'-lookup_table',$lut;
    push @args,'-discrete' if $discrete;
    
  } else {
    push @args,'-spectral';
  }
  
  push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
  do_cmd(@args);
  
  #mix colors
  if($over) {
    if($#ovl_range>0)
    {
      do_cmd('minccalc','-express',"A[0]>$ovl_range[0]?1:0",'-byte',"$tmpdir/overlay.mnc","$tmpdir/ovl_mask.mnc");

      do_cmd('minclookup','-gray',"$tmpdir/ovl_mask.mnc","$tmpdir/ovl_mask_g.mnc",'-min',0,'-max',1);

      do_cmd('minccalc','-express','A[2]>0.5?A[1]:A[0]',
           "$tmpdir/gray.mnc","$tmpdir/overlay_rgb.mnc","$tmpdir/ovl_mask_g.mnc","$tmpdir/gray_.mnc");
      
    } else {
      print STDERR "Use --ovl-range with --over!\n";
    }
    
  } elsif($max) {
    do_cmd('mincmath','-max',"$tmpdir/gray.mnc","$tmpdir/overlay_rgb.mnc","$tmpdir/gray_.mnc");
  } else {
    do_cmd('mincaverage',"$tmpdir/gray.mnc","$tmpdir/overlay_rgb.mnc","$tmpdir/gray_.mnc");
  }
  do_cmd('mv',"$tmpdir/gray_.mnc","$tmpdir/gray.mnc");
}

if($red || $green || $blue)
{
  my @channels;
  my @channels_masks;
  if($red)
  {
    do_cmd('mincresample',$red,'-like',$in,"$tmpdir/red.mnc");
    if($mask)
    {
      do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/red.mnc",$mask,"$tmpdir/masked_red.mnc");

      do_cmd('minccalc','-express',"A[0]>0.5&&A[1]>$ovl_range[0]?1:0",'-byte',$mask,"$tmpdir/red.mnc","$tmpdir/ovl_mask_red.mnc") if $over && $#ovl_range>0;
      
      do_cmd('mv',"$tmpdir/masked_red.mnc","$tmpdir/red.mnc");
    } else {
      do_cmd('minccalc','-express',"A[0]>$ovl_range[0]?1:0",'-byte',"$tmpdir/red.mnc","$tmpdir/ovl_mask_red.mnc") if $over && $#ovl_range>0;
    }
    @args=("minclookup","$tmpdir/red.mnc","$tmpdir/overlay_red.mnc",'-byte','-lut_string','0.0 0.0 0.0 0.0;1.0 1.0 0.0 0.0');
    push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
    do_cmd(@args);
    push @channels,"$tmpdir/overlay_red.mnc";
    push @channels_masks,"$tmpdir/ovl_mask_red.mnc" if $over && $#ovl_range>0;
  }
  
  if($green)
  {
    do_cmd('mincresample',$green,'-like',$in,"$tmpdir/green.mnc");
    if($mask)
    {
      do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/green.mnc",$mask,"$tmpdir/masked_green.mnc");

      do_cmd('minccalc','-express',"A[0]>0.5&&A[1]>$ovl_range[0]?1:0",'-byte',$mask,"$tmpdir/green.mnc","$tmpdir/ovl_mask_green.mnc") if $over && $#ovl_range>0;
      
      do_cmd('mv',"$tmpdir/masked_green.mnc","$tmpdir/green.mnc");
    } else {

      do_cmd('minccalc','-express',"A[0]>$ovl_range[0]?1:0",'-byte',"$tmpdir/green.mnc","$tmpdir/ovl_mask_green.mnc") if $over && $#ovl_range>0;
    }
    @args=("minclookup","$tmpdir/green.mnc","$tmpdir/overlay_green.mnc",'-byte','-lut_string','0.0 0.0 0.0 0.0;1.0 0.0 1.0 0.0');
    push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
    do_cmd(@args);
    push @channels,"$tmpdir/overlay_green.mnc";
    push @channels_masks,"$tmpdir/ovl_mask_green.mnc" if $over && $#ovl_range>0;
  }
  
  if($blue)
  {
    do_cmd('mincresample',$blue,'-like',$in,"$tmpdir/blue.mnc");
    if($mask)
    {
      do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',"$tmpdir/blue.mnc",$mask,"$tmpdir/masked_blue.mnc");
      do_cmd('minccalc','-express',"A[0]>0.5&&A[1]>$ovl_range[0]?1:0",'-byte',$mask,"$tmpdir/blue.mnc","$tmpdir/ovl_mask_blue.mnc") if $over && $#ovl_range>0;
      
      do_cmd('mv',"$tmpdir/masked_blue.mnc","$tmpdir/blue.mnc");
    } else {
      do_cmd('minccalc','-express',"A[0]>$ovl_range[0]?1:0",'-byte',"$tmpdir/blue.mnc","$tmpdir/ovl_mask_blue.mnc") if $over && $#ovl_range>0;
    }
    @args=("minclookup","$tmpdir/blue.mnc","$tmpdir/overlay_blue.mnc",'-byte','-lut_string','0.0 0.0 0.0 0.0;1.0 0.0 0.0 1.0');
    push @args,'-min',$ovl_range[0],'-max',$ovl_range[1] if $#ovl_range>0;
    do_cmd(@args);
    push @channels,"$tmpdir/overlay_blue.mnc";
    push @channels_masks,"$tmpdir/ovl_mask_blue.mnc" if $over && $#ovl_range>0;
  }
   
  
  #mix colors
  if($over && $#channels_masks>-1 ) {
    do_cmd('mincmath','-max',@channels_masks,"$tmpdir/ovl_mask.mnc");
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
@args=('mincpik',"$tmpdir/gray.mnc","$tmpdir/coronal.miff",'-coronal','-slice',$coronal);
push @args,'-scale',$scale if $scale;
do_cmd(@args);

@args=('mincpik',"$tmpdir/gray.mnc","$tmpdir/axial.miff",'-axial','-slice',$axial);
push @args,'-scale',$scale if $scale;
do_cmd(@args);

@args=('mincpik',"$tmpdir/gray.mnc","$tmpdir/sagittal.miff",'-sagittal','-slice',$sagittal);
push @args,'-scale',$scale if $scale;
do_cmd(@args);

if($separate)
{
  do_cmd("convert","$tmpdir/coronal.miff",$out_coronal);
  do_cmd("convert","$tmpdir/axial.miff",$out_axial);
  do_cmd("convert","$tmpdir/sagittal.miff",$out_sagittal);
} else {
  my $geo='+0+0';

  if($vertical)
  {
    equalize_width("$tmpdir/axial.miff",
         "$tmpdir/sagittal.miff",
         "$tmpdir/coronal.miff");
  } else {
    equalize_height("$tmpdir/axial.miff",
         "$tmpdir/sagittal.miff",
         "$tmpdir/coronal.miff");
  }
  #chomp($geo);
  my @args=('montage','-depth',8,#'-texture','rose:',
          '-tile',$vertical?'1x3':'3x1',
          '-geometry',$geo,'-gravity','North',
          '-background' ,$background,
          '-bordercolor',$background,
          '-mattecolor',$background,
          '-fill',$foreground,
          '-stroke',$foreground,
         "$tmpdir/axial.miff",
         "$tmpdir/sagittal.miff",
         "$tmpdir/coronal.miff");
  push @args,'-pointsize',$pointsize if $pointsize;
  push @args,'-title',$title if $title;
  push @args,$out;
  
  do_cmd(@args);
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
