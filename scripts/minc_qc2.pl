#!/usr/bin/env perl

##############################################################################
# 
# pipeline_qc.pl
#
# Input:
#
# Output:
#      o a jpeg file with the mask overlaid on the T1 image.
#      o a jpeg file with the T1 image overlaid on the T2 image.
#
#
# Larry Baer, May, 2005
# McConnell Brain Imaging Centre, 
# Montreal Neurological Institute, 
# McGill University
##############################################################################

use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;

my $me = basename($0);
my $verbose = 0;
my $fake = 0;
my $clobber = 0;
my $title;
my $mask;
my $spectral;
my @image_range;
my $spectral_mask;
my @mask_range;
my $horizontal=0;
my $big=0;
my $avg=0;
my $bg='grey10';

GetOptions(
	   'verbose' => \$verbose,
	   'fake'    => \$fake,
	   'clobber' => \$clobber,
	   'title=s'          => \$title,
	   'mask=s' 	      => \$mask,
	   'spectral'         => \$spectral,
	   'spectral-mask'    => \$spectral_mask,
	   'image-range=f{2}' => \@image_range,
	   'mask-range=f{2}' => \@mask_range,
     'horizontal'      => \$horizontal,
     'avg'             => \$avg,
     'big'             => \$big,
     'bg=s'            => \$bg
	   );

die "Usage: $me  <input.mnc> <output jpeg file> [--title <title> --clobber --verbose --mask <file> --spectral --image-range a b --spectral-mask --mask-range a b]\n"  if $#ARGV < 1;

#########################
# Get the arguments.
my $infile;

my $out=pop @ARGV;
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );


delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};

#print $#image_range,"\n";

foreach $infile(@ARGV)
{
  # make tmpdir

  my $outfile;
  unless($#ARGV) #we should have directory
  {
    $outfile=$out;
  } else {
    $outfile=$out.'/'.basename($infile).'.jpg';
    $title=basename($infile);
  }
  if($clobber || check_file($outfile) )
  {
    my @args=('minclookup','-clobb',$spectral?'-spectral':'-grey',#'-range',10,100,
      $infile, "$tmpdir/grey.mnc",'-clobber');

    if($#image_range>0 && $image_range[0]!=$image_range[1])
    {
      push @args,'-range',$image_range[0],$image_range[1];
    }
    
    do_cmd(@args);
    
    #red mask
    if($mask)
    {
      my @range=split(/\n/,`mincstats -q -min -max $mask`);
      my @args=('minclookup','-clobb',$mask, "$tmpdir/mask.mnc",'-clobber','-range',$range[0],$range[1]);
      if($spectral_mask) {
        push @args,'-spectral'
      } else {
        push @args,'-lut_string','0.0 0.0 0.0 0.0;1.0 1.0 0.0 0.0';
      }
      if($#mask_range>0 && $mask_range[0]!=$mask_range[1])
      {
        push @args,'-range',$mask_range[0],$mask_range[1];
      }
      do_cmd(@args);
      
      unless($avg)
      {
        do_cmd('mincmath','-max',"$tmpdir/grey.mnc","$tmpdir/mask.mnc","$tmpdir/masked.mnc");
      } else {
        do_cmd('mincaverage',"$tmpdir/grey.mnc","$tmpdir/mask.mnc","$tmpdir/masked.mnc");
      }
      do_cmd('mv',"$tmpdir/masked.mnc","$tmpdir/grey.mnc");
    }
  
    make_multipane("$tmpdir/grey.mnc", $title, $outfile);
  }
}


#####################################################################
# sub-routine to make a multipane view 
sub make_multipane{
   my($mncfile, $text, $imgfile, @ext_args) = @_;
   my $smalltilesize = $big?250:150;
   my(@args, @mont_args,@mont_x,@mont_y,@mont_z);
   
   # try a .gz if file missing
   $mncfile .= '.gz' if (!-e $mncfile);
   my($xspace,$yspace,$zspace)=split(/\n/,`mincinfo -dimlength xspace  -dimlength yspace  -dimlength zspace $mncfile`);
   $xspace=$xspace*1.0;
   $yspace=$yspace*1.0;
   $zspace=$zspace*1.0;
   #z 
   foreach  (30,53,76,99,122,145) {
       @args = ('mincpik', '-scale','1','-transverse','-slice',int($zspace*$_/189.0),'-clobber'); #linux:add -clobber
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/T$_.miff");
       &do_cmd(@args);
       
       push(@mont_z, "$tmpdir/T$_.miff");
       }
       
  &do_cmd('montage',
          '-tile', $horizontal?'2x3':'3x2',
          '-background', $bg,
          '-geometry', '+0+0',
          @mont_z,
          "$tmpdir/mont_z.miff");
          
  push @mont_args,"$tmpdir/mont_z.miff";
 #x
 foreach  (50,60,70,130,120,110) {
       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$_/197.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$_.miff");
       &do_cmd(@args);
       
       push(@mont_x, "$tmpdir/S$_.miff");
       }

  &do_cmd('montage',
          '-tile', $horizontal?'2x3':'3x2',
          '-background', $bg,
          '-geometry', '+0+0',
          @mont_x,
          "$tmpdir/mont_x.miff");
  push @mont_args,"$tmpdir/mont_x.miff";
 #y
 foreach  (60,80,110,120,140,160) {
     @args = ('mincpik','-scale','1','-coronal','-slice',int($xspace*$_/233.0),'-clobber');
     push(@args, @ext_args) if @ext_args;
     push(@args, $mncfile, "$tmpdir/C$_.miff");
     &do_cmd(@args);
     
     push(@mont_y, "$tmpdir/C$_.miff");
     }
     
  &do_cmd('montage',
          '-tile', $horizontal?'2x3':'3x2',
          '-background', $bg,
          '-geometry', '+0+0',
          @mont_y,
          "$tmpdir/mont_y.miff");
  push @mont_args,"$tmpdir/mont_y.miff";
 
    # do the montage
  &do_cmd('montage',
          '-tile', ($horizontal?'3x1':'1x3'),
          '-background', $bg,
          '-gravity','center',
          '-geometry', '+0+0',
          @mont_args,
          "$tmpdir/mont.miff");
             
    # Add the title
  my $w=$horizontal?$smalltilesize*6:$smalltilesize*3;
  my $h=$horizontal?$smalltilesize*3:$smalltilesize*6;
  
  @args = ('convert', '-box', 'white',
               '-resize',"${w}x${h}");
     #'-font', '7x13', 
     #'-fill', 'white',
     #'-draw', "text 2,15 \"$text\"")
   push(@args,'-draw', "text 2,15 \"$text\"") if $text;
   #push(@args, @more_args) if @more_args;
  &do_cmd(@args,"$tmpdir/mont.miff", $imgfile);
    
 }



#####################################################################
sub do_cmd { 
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
      system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  if(-e $_[0])
  {
    warn("${_[0]} exists!\n");
    return 0;
  }    
  return 1;
}