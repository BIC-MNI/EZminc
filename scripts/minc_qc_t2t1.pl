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
my @t1_range=(10,60);
my @t2_range=(70,90);

GetOptions(
	   'verbose' => \$verbose,
	   'fake'    => \$fake,
	   'clobber' => \$clobber,
     'title=s' => \$title,
     'mask=s'  => \$mask
	   );

die "Usage: $me <input_t1.mnc> <input_t2.mnc> <output jpeg file> [--title <title> --clobber --verbose ]\n"  if $#ARGV < 1;

#########################
# Get the arguments.
my $infile;

my $out=pop @ARGV;
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );


delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};

my ($in_t1,$in_t2)=@ARGV;

{
  # make tmpdir

  my $outfile;
  $outfile=$out;
  if($clobber || check_file($outfile) )
  {
    
    do_cmd('minclookup','-clobb','-lut_string','0.0 0.0 0.0 0.0;1.0 1.0 0.0 0.0','-range',$t1_range[0],$t1_range[1],
        $in_t1, "$tmpdir/red.mnc",'-clobber');
    
    do_cmd('minclookup','-clobb','-lut_string','0.0 0.0 0.0 0.0;1.0 0.0 1.0 0.0','-range',$t2_range[0],$t2_range[1],
        $in_t2, "$tmpdir/green.mnc",'-clobber');
    
    do_cmd('mincmath','-max',"$tmpdir/red.mnc","$tmpdir/green.mnc","$tmpdir/comp.mnc");
 
    make_multipane("$tmpdir/comp.mnc", $title, $outfile);
  }
}


#####################################################################
# sub-routine to make a multipane view 
sub make_multipane{
   my($mncfile, $text, $imgfile, @ext_args) = @_;
   my $smalltilesize = 150;
   my(@args, @mont_args);
   
   # try a .gz if file missing
   $mncfile .= '.gz' if (!-e $mncfile);
   my($xspace,$yspace,$zspace)=split(/\n/,`mincinfo -dimlength xspace  -dimlength yspace  -dimlength zspace $mncfile`);
   $xspace=$xspace*1.0;
   $yspace=$yspace*1.0;
   $zspace=$zspace*1.0;
   #z 
   foreach  (30.0,35.0,40.0,45.0,50.0,145.0) {
       @args = ('mincpik', '-scale','1','-transverse','-slice',int($zspace*$_/189.0),'-clobber'); #linux:add -clobber
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/T$_.miff");
       &do_cmd(@args);
       
       push(@mont_args, "$tmpdir/T$_.miff");
       }
 #x
 foreach  (50,60,70,130,120,110) {
       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$_/197.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$_.miff");
       &do_cmd(@args);
       
       push(@mont_args, "$tmpdir/S$_.miff");
       }
 #y
 foreach  (60,80,110,120,140,160) {
     @args = ('mincpik','-scale','1','-coronal','-slice',int($xspace*$_/233.0),'-clobber');
     push(@args, @ext_args) if @ext_args;
     push(@args, $mncfile, "$tmpdir/C$_.miff");
     &do_cmd(@args);
     
     push(@mont_args, "$tmpdir/C$_.miff");
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