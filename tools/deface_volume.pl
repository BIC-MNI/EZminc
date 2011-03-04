#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $amp=12;
my $fwhm=6;
my $keep_tmp=0;
my $tal_grid;
#my $mask_distance=5;
my $edge_smooth=2;
my $tal_xfm;
my $save_grid;
my $model;
my $model_dir;
my $face;
my $brain;
my $normalize;

GetOptions (    
          "verbose"       => \$verbose,
          "clobber"       => \$clobber,
          "fwhm=f"        => \$fwhm,
          "amp=f"         => \$amp,
          "tal_grid=s"    => \$tal_grid,  
          "keep_tmp"      => \$keep_tmp,
          "edge_smooth=f" => \$edge_smooth,
          "model=s"       => \$model,
          "model_dir=s"   => \$model_dir,
          "tal_xfm=s"     => \$tal_xfm,
          "save_grid=s"   => \$save_grid,
          "brain=s"       => \$brain,
          "face=s"        => \$face,
          "normalize"     => \$normalize,
          ); 

die <<HELP
Usage: $me <scan_in> <output> 
 [ --verbose 
   --clobber 
   --fwhm <f> 
   --amp <f> 
   --tal_grid <grid> 
   --keep_tmp 
   --edge_smooth <f> 
   --save_grid <random_grid> 
   --tal_xfm <xfm> 
   --model <name> 
   --model_dir <dir> 
   --face <face> 
   --brain <brain>
 ]
HELP
if $#ARGV<1;

my ($scan,$output)=@ARGV;

check_file($output) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => !$keep_tmp );

#if not registered - register first
#assume it is not normalized
if($model && $model_dir && !$tal_xfm)
{
  my $modelfn="$model_dir/${model}.mnc";
  my $in=$scan;
  if($normalize)
  {
    do_cmd('nu_correct'  , "-clobber", "-iter", 100, "-stop", 0.0001, "-fwhm", 0.1,$scan, "$tmpdir/nuc.mnc");
    do_cmd('volume_pol'  , '--order', 1, '--min', 0, '--max', 100,  "$tmpdir/nuc.mnc", $modelfn,'--expfile', "$tmpdir/stats", '--clobber');
    do_cmd('minccalc'    , "$tmpdir/nuc.mnc", "$tmpdir/clamp.mnc", '-expfile', "$tmpdir/stats", '-clobber','-short');
    $in="$tmpdir/clamp.mnc";
  }
  do_cmd('bestlinreg_s','-lsq9',$in,$modelfn,"$tmpdir/tal.xfm");
  $tal_xfm="$tmpdir/tal.xfm";
}

if($face && $brain)
{
  do_cmd('itk_morph','--exp','D[2]',$brain,"$tmpdir/brain.mnc");
  do_cmd('itk_morph','--exp','D[2]',$face ,"$tmpdir/face.mnc");
  #substract brain from face
  do_cmd('minccalc','-expression','A[0]==1&&A[1]==0?1:0', "$tmpdir/face.mnc", "$tmpdir/brain.mnc", "$tmpdir/mask.mnc",'-clobber');
  do_cmd('itk_morph','--exp','E[1]',"$tmpdir/mask.mnc","$tmpdir/face.mnc",'--clobber');
  $face="$tmpdir/face.mnc";
} else {
  $face="$model_dir/${model}_face_mask.mnc" unless $face;
}

#create random grid
if(!$tal_grid)
{
  #dilate brain,face
  do_cmd('make_random_grid.pl', $face, '--edge_smooth', $edge_smooth, '--mask', $face, "$tmpdir/deface_grid.mnc",'--amplitude',$amp,'--fwhm',$fwhm);
  $tal_grid="$tmpdir/deface_grid.mnc";
  do_cmd('cp',$tal_grid,$save_grid) if $save_grid;
}

#transform to native space
do_cmd('mincresample', '-nearest', $face, '-transform', $tal_xfm, '-invert_transformation', "$tmpdir/native_face.mnc",  '-like', $scan);

do_cmd('xfminvert',$tal_xfm,"$tmpdir/tal_to_native.xfm");
do_cmd('itk_resample',$scan,"$tmpdir/scan.mnc",'--uniformize',2,'--order',0);
do_cmd('resample_grid',$tal_grid,"$tmpdir/tal_to_native.xfm","$tmpdir/random_grid_0.mnc",'--like',"$tmpdir/scan.mnc");

open XFM,">$tmpdir/random.xfm" or die;
print XFM "MNI Transform File\nTransform_Type = Grid_Transform;\nDisplacement_Volume = random_grid_0.mnc;\n";
close XFM;
#modulate scan
do_cmd('mincresample', '-nearest', '-transform', "$tmpdir/random.xfm", $scan, $output, '-clobber', '-use_input_sampling');

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}
sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
