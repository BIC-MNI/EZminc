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
my $cyanred;
my $cyanred_mask;
my $lut;
my $gray_mask;
my $discrete;
my $hotmetal_mask;
my $big;
my $bbox;
my $labels_lut;
my $labels;
my $labels_mask;
my $red;
my $green_mask;
my $clamp;
my $mask_lut;
my $discrete_mask;
my $slices=1;

GetOptions(
      'verbose' => \$verbose,
      'fake'    => \$fake,
      'clobber' => \$clobber,
      'slices=n' => \$slices
      );

die "Usage: $me  <input.mnc> <output prefix> 
[ 
  --clobber 
  --verbose 
  --slices <n>
]\n"  if $#ARGV < 1;

#########################
# Get the arguments.

my $out_prefix=pop @ARGV;
my $infile=pop @ARGV;
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};

  # make tmpdir
  
my @args=('minclookup','-clobb', $infile, "$tmpdir/grey.mnc",'-clobber','-gray');

do_cmd(@args);

make_multipane_big("$tmpdir/grey.mnc", $title, $out_prefix,$slices);



#####################################################################
# sub-routine to make a multipane view 
sub make_multipane_big{
    my($mncfile, $text, $imgfile, $slices, @ext_args) = @_;
    my $smalltilesize = 170;
    my $i;
    my (@args, @piks);

    my $steps=20;
    my $columns=$steps/2;
    my $rows=6;
    my ($x,$y,$z);
    # try a .gz if file missing
    $mncfile .= '.gz' if (!-e $mncfile);
    
    my($xspace,$yspace,$zspace)=split(/\n/,`mincinfo -dimlength xspace -dimlength yspace -dimlength zspace $mncfile`);
    $xspace=$xspace*1.0;
    $yspace=$yspace*1.0;
    $zspace=$zspace*1.0;
    
    my @zslices=(1.0/2);
    my @yslices=(1.0/2);
    my @xslices=(1.0/2);
    # HACK
    if($slices>1) {
        @zslices=(53.0/192, 111.0/192.0);
    }
    #z
    foreach $z(@zslices)  {
        @args = ('mincpik', '-scale', '1', '-transverse', '-slice',int($zspace*$z),'-clobber','-verbose'); 
        push(@args, @ext_args) if @ext_args;
        push(@args, $mncfile, "$tmpdir/T_$z.miff" );
        do_cmd(@args);
        push(@piks,"$tmpdir/T_$z.miff");
    }
    #x
    foreach $x (@xslices) {
        @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$x),'-clobber','-verbose');
        push(@args, @ext_args) if @ext_args;
        push(@args, $mncfile, "$tmpdir/S_$x.miff");
        do_cmd(@args);
        push(@piks,"$tmpdir/S_$x.miff");
    }
    #y
    foreach $y (@yslices) {
        @args = ('mincpik','-scale','1','-coronal','-slice',int($yspace*$y),'-clobber','-verbose');
        push(@args, @ext_args) if @ext_args;
        push(@args, $mncfile, "$tmpdir/C_$y.miff");
        do_cmd(@args);
        push(@piks,"$tmpdir/C_$y.miff");
    }    
    
    for($i=0;$i<=$#piks;$i=$i+1)
    {
        do_cmd('convert','-background', 'black', '-gravity', 'center', '-resize', '256x256',
                         '-gravity','center', '-extent', '224x224', $piks[$i], "${imgfile}_${i}.jpg");
    }
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

