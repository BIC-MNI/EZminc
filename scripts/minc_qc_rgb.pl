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
my @image_range;
my $big;

GetOptions(
	   'verbose' => \$verbose,
	   'fake'    => \$fake,
	   'clobber' => \$clobber,
     'title=s' => \$title,
     'mask=s'  => \$mask,
     'big'             => \$big,
	   'image-range=f{2}' => \@image_range,
	   );

die "Usage: $me <input_r.mnc> [input_g.mnc] [input_b.mnc] <output jpeg file> [--title <title> --clobber --verbose --image-range a b ]\n"  if $#ARGV < 1;

#########################
# Get the arguments.
my $infile;

my $out=pop @ARGV;
check_file($out) unless $clobber;

# make tmpdir
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};
my @luts=('0.0 0.0 0.0 0.0;1.0 1.0 0.0 0.0', '0.0 0.0 0.0 0.0;1.0 0.0 1.0 0.0', '0.0 0.0 0.0 0.0;1.0 0.0 0.0 1.0');

my $i;

for($i=0;$i<=$#ARGV;$i+=1 )
{
  my @args=('minclookup','-clobb','-lut_string',$luts[$i],#'-range',$t1_range[0],$t1_range[1],
      $ARGV[$i], "$tmpdir/cur.mnc",'-clobber');
      
  if($#image_range>0 && $image_range[0]!=$image_range[1])
    {
      push @args,'-range',$image_range[0],$image_range[1];
    }
    
  do_cmd(@args);
  
  if( -e "$tmpdir/comp.mnc") {
    do_cmd('mincmath','-max',"$tmpdir/cur.mnc","$tmpdir/comp.mnc","$tmpdir/comp_.mnc");
    do_cmd('mv',"$tmpdir/comp_.mnc","$tmpdir/comp.mnc")
  } else {
    do_cmd('mv',"$tmpdir/cur.mnc","$tmpdir/comp.mnc")
  }
}

if($big) {
  make_multipane_big("$tmpdir/comp.mnc", $title, $out);
} else {
  make_multipane("$tmpdir/comp.mnc", $title, $out);
}

#####################################################################
# sub-routine to make a multipane view 
sub make_multipane{
   my($mncfile, $text, $imgfile, @ext_args) = @_;
   my $smalltilesize = 150;
   my(@args, @mont_args);
   
   # try a .gz if file missing
   $mncfile .= '.gz' if (!-e $mncfile);
   my($xspace,$yspace,$zspace)=split(/\n/,`mincinfo -dimlength xspace -dimlength yspace  -dimlength zspace $mncfile`);
   $xspace=$xspace*1.0;
   $yspace=$yspace*1.0;
   $zspace=$zspace*1.0;
   my $i;
   #z 
   foreach  $i(30.0,35.0,40.0,45.0,50.0,145.0) {
       @args = ('mincpik', '-scale','1','-transverse','-slice',int($zspace*$i/181.0),'-clobber'); #linux:add -clobber
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/T${i}.miff" );
			 do_cmd(@args);
       
       push(@mont_args, "$tmpdir/T${i}.miff");
       }
 #x
  foreach  $i(50,60,70,130,120,110) {
       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$i/181.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$i.miff");
       do_cmd(@args);
       
       push(@mont_args, "$tmpdir/S$i.miff");
       }
 #y
  foreach  $i(60,80,110,120,140,160) {
     @args = ('mincpik','-scale','1','-coronal','-slice',int($yspace*$i/217.0),'-clobber');
     push(@args, @ext_args) if @ext_args;
     push(@args, $mncfile, "$tmpdir/C$i.miff");
     &do_cmd(@args);
     
     push(@mont_args, "$tmpdir/C$i.miff");
     }
 
    # do the montage
  &do_cmd('montage',
          '-tile', '3x6',
          '-background', 'grey10',
          '-geometry', $smalltilesize . 'x' . $smalltilesize . '+1+1',
          @mont_args,
          "$tmpdir/mont.miff");
             
    # Add the title
  @args = ('convert', '-box', 'white');
     #'-font', '7x13', 
     #'-fill', 'white',
     #'-draw', "text 2,15 \"$text\"");
   push(@args,'-draw', "text 2,15 \"$text\"") if $text;
   #push(@args, @more_args) if @more_args;
  &do_cmd(@args,"$tmpdir/mont.miff", $imgfile);
    
 }

#####################################################################
# sub-routine to make a multipane view 
sub make_multipane_big{
   my($mncfile, $text, $imgfile, @ext_args) = @_;
   my $smalltilesize = 170;
   my(@args, @mont_args);

   my $steps=20;
   my $columns=$steps/2;
   my $rows=6;

   # try a .gz if file missing
   $mncfile .= '.gz' if (!-e $mncfile);
   my($xspace,$yspace,$zspace)=split(/\n/,`mincinfo -dimlength xspace -dimlength yspace  -dimlength zspace $mncfile`);
   $xspace=$xspace*1.0;
   $yspace=$yspace*1.0;
   $zspace=$zspace*1.0;
   my $i;
   my $j;
   #z 
   for($j=0;$j<$steps;$j++) {

       $i=int(10+(150.0-10.0)*$j/($steps-1));

       @args = ('mincpik', '-scale','1','-transverse','-slice',int($zspace*$i/181.0),'-clobber'); #linux:add -clobber
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/T${i}.miff" );
			 do_cmd(@args);
       
       push(@mont_args, "$tmpdir/T${i}.miff");
       }
 #x
   for($j=0;$j<($steps/2);$j++) {

       $i=int(28.0+(166.0-28.0)*$j/($steps-1));

       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$i/193.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$i.miff");
       do_cmd(@args);
       
       push(@mont_args, "$tmpdir/S$i.miff");
       }

   for($j=($steps-1);$j>=($steps/2);$j--) {

       $i=int(28.0+(166.0-28.0)*$j/($steps-1));

       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$i/193.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$i.miff");
       do_cmd(@args);
       
       push(@mont_args, "$tmpdir/S$i.miff");
       }
 #y
   for($j=0;$j<$steps;$j++) {

     $i=int(25+(195.0-25.0)*$j/($steps-1));
     @args = ('mincpik','-scale','1','-coronal','-slice',int($yspace*$i/217.0),'-clobber');
     push(@args, @ext_args) if @ext_args;
     push(@args, $mncfile, "$tmpdir/C$i.miff");
     &do_cmd(@args);
     
     push(@mont_args, "$tmpdir/C$i.miff");
     }
 
    # do the montage
  &do_cmd('montage',
          '-tile', "${columns}x${rows}",
          '-background', 'grey10',
          '-geometry', $smalltilesize . 'x' . $smalltilesize . '+1+1',
          @mont_args,
          "$tmpdir/mont.miff");
             
    # Add the title
  @args = ('convert', '-box', 'white');
     #'-font', '7x13', 
     #'-fill', 'white',
     #'-draw', "text 2,15 \"$text\"");
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
    die("${_[0]} exists!\n");
    return 0;
  }    
  return 1;
}