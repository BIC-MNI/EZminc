#! /usr/bin/env perl


############################# MNI Header #####################################
#@NAME       :  phantomfit.pl
#@DESCRIPTION:  estimate distortion field based on nonlinear registration using ANTs
#@COPYRIGHT  :
#              Vladimir S. Fonov  February, 2012
#              Montreal Neurological Institute, McGill University.
#              Permission to use, copy, modify, and distribute this
#              software and its documentation for any purpose and without
#              fee is hereby granted, provided that the above copyright
#              notice appear in all copies.  The author and McGill University
#              make no representations about the suitability of this
#              software for any purpose.  It is provided "as is" without
#              express or implied warranty.
###############################################################################


use strict;
use warnings "all";
use Getopt::Tabular;
use File::Basename;
use File::Temp qw/ tempdir /;


my $elastix_fast= <<ELX1;
(FixedInternalImagePixelType "float")
(MovingInternalImagePixelType "float")
(FixedImageDimension 3)
(MovingImageDimension 3)
(UseDirectionCosines "true")

(Registration "MultiResolutionRegistration")
(Interpolator "BSplineInterpolator" )
(ResampleInterpolator "FinalBSplineInterpolator" )
(Resampler "DefaultResampler" )

(FixedImagePyramid  "FixedSmoothingImagePyramid")
(MovingImagePyramid "MovingSmoothingImagePyramid")

(Optimizer "AdaptiveStochasticGradientDescent")
(Transform "BSplineTransform")
(Metric "AdvancedNormalizedCorrelation")

(FinalGridSpacingInPhysicalUnits 12)

(HowToCombineTransforms "Compose")

(ErodeMask "false")

(NumberOfResolutions 3)

(ImagePyramidSchedule 8 8 8  4 4 4  2 2 2)

(MaximumNumberOfIterations 1000)
(MaximumNumberOfSamplingAttempts 3)

(NumberOfSpatialSamples 4096)

(NewSamplesEveryIteration "true")
(ImageSampler "Random" )

(BSplineInterpolationOrder 1)

(FinalBSplineInterpolationOrder 3)

(DefaultPixelValue 0)

(WriteResultImage "false")

// The pixel type and format of the resulting deformed moving image
(ResultImagePixelType "float")
(ResultImageFormat "mnc")
ELX1

my $elastix_par= <<ELX;
(FixedInternalImagePixelType "float")
(MovingInternalImagePixelType "float")
(FixedImageDimension 3)
(MovingImageDimension 3)
(UseDirectionCosines "true")

(Registration "MultiResolutionRegistration")
(Interpolator "BSplineInterpolator" )
(ResampleInterpolator "FinalBSplineInterpolator" )
(Resampler "DefaultResampler" )

(FixedImagePyramid  "FixedSmoothingImagePyramid")
(MovingImagePyramid "MovingSmoothingImagePyramid")

(Optimizer "AdaptiveStochasticGradientDescent")
(Transform "BSplineTransform")
(Metric "AdvancedNormalizedCorrelation")

(FinalGridSpacingInPhysicalUnits 12)

(HowToCombineTransforms "Compose")

(ErodeMask "false")

(NumberOfResolutions 4)

(ImagePyramidSchedule 8 8 8  4 4 4  2 2 2 1 1 1 )

(MaximumNumberOfIterations 2000 2000 2000 8000)
(MaximumNumberOfSamplingAttempts 3)

(NumberOfSpatialSamples 4096)

(NewSamplesEveryIteration "true")
(ImageSampler "Random" )

(BSplineInterpolationOrder 1)

(FinalBSplineInterpolationOrder 3)

(DefaultPixelValue 0)

(WriteResultImage "false")

// The pixel type and format of the resulting deformed moving image
(ResultImagePixelType "float")
(ResultImageFormat "mnc")
ELX

my($Help, $Usage, $me);
my(@opt_table, %opt, $outxfm, $outfile, @args, $tmpdir);

$me = &basename($0);
%opt = (
   'verbose'   => 1,
   'debug'     => 0,
   'clobber'   => 0,
   'fake'      => 0,
   'init_xfm'  => undef,
   'order'     => 5,
   'par'       => undef,
   'measure'   => undef,
   'step_iterations'=> undef,
   'min_step'  => 2,
   'limit'     => 0,
   'weight'    => 1,
   'keep'      => 1.0,
   'cyl'       => 0,
   'init'      => undef,
   'work_dir'  => undef,
   'pca'       => undef,
   'pcs'       => undef,
   );

$Help = <<HELP;
| $me does non-linear fitting between N files
| restricted by spherical harmonics functions
| 
| Problems or comments should be sent to: vladimir.fonov\@gmail.com
HELP

$Usage = "Usage: $me [options] source.mnc target.mnc source_mask.mnc target_mask.mnc bricks.mnc [source.mnc target.mnc source_mask.mnc target_mask.mnc bricks.mnc...] output.xfm\n$me -help to list options\n\n";

@opt_table = (
   ["-verbose", "boolean", 0, \$opt{verbose},
      "be verbose" ],
   ["-debug", "boolean", 0, \$opt{debug},
      "for debugging" ],
   ["-clobber", "boolean", 0, \$opt{clobber},
      "clobber existing check files" ],
   ["-fake", "boolean", 0, \$opt{fake},
      "do a dry run, (echo cmds only)" ],
   ["-init_xfm", "string", 1, \$opt{init_xfm},
      "initial transformation (default identity)" ],
   ["-order","integer",1,\$opt{order},
      "Spherical harmonics order"],
   ["-par", "string", 1, \$opt{par},
      "Output parameters into a file" ],
   ["-measure", "string", 1, \$opt{measure},
      "Output measurements into a file" ], 
   ["-step_iterations", "integer", 1, \$opt{step_iterations},
      "Maximum number of iterations per step" ], 
   ["-min_step", "float", 1, \$opt{min_step},
      "Minimal step size for nonlinear registration (min 1.0)" ], 
   ["-limit", "boolean", 0, \$opt{limit},
      "Limit linear component to identity" ],
   ["-keep", "float", 1, \$opt{keep},
      "Fraction of points to keep for LTSQ algorithm (0-1]" ],
   ["-cylindric", "boolean", 0, \$opt{cyl},
      "Use cylindric functions instead of spherical" ],
   ["-init", "string", 1, \$opt{init},
      "Initial estimation" ],
   ["-work_dir", "string", 1, \$opt{work_dir},
      "Work directory (instead of temp), usefull for debugging" ],
   ["-pca", "string", 1, \$opt{pca},
      "Use pca rotation matrix" ],
   ["-pcs", "integer", 1, \$opt{pcs},
      "limit number of PCs" ],
   ["-outdir","string", 1, \$opt{outdir},
      "Output files directory (bricks)" ],
   );

# Check arguments
&Getopt::Tabular::SetHelp($Help, $Usage);
&GetOptions (\@opt_table, \@ARGV) || exit 1;
die $Usage if  $#ARGV < 4;
my @source;
my @target;
my @source_mask;
my @target_mask;
my @bricks;

print "ARGV:",join(',',@ARGV),"\n";

my $minc_compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $minc_compress;

$outxfm = pop(@ARGV);
$outfile = $opt{output};

for(my $i=0;$i <= (($#ARGV-1)/5);$i++)
{
  push @source,$ARGV[$i*5];
  push @target,$ARGV[$i*5+1];
  push @source_mask,$ARGV[$i*5+2];
  push @target_mask,$ARGV[$i*5+3];
  push @bricks,$ARGV[$i*5+4];
}
print "Sources:",join(' ',@source),"\n";
print "Target:",join(' ',@target),"\n";

check_file($outxfm)   unless $opt{clobber};
check_file($outfile)  unless $opt{clobber} || !defined($outfile);
check_file($opt{par}) unless $opt{clobber} || !defined($opt{par});

# make tmpdir
unless($opt{work_dir})
{
  $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
} else {
  $tmpdir = $opt{work_dir};
  do_cmd('mkdir','-p',$tmpdir);
}

# set up filename base
my($i, $s_base, $t_base, $tmp_xfm, $tmp_grid, $tmp_source, $tmp_target, $prev_grid, $prev_xfm);
my @grids;
my @masks;
# a fitting we shall go...
# TODO: make this parallel?
for(my $k=0;$k<=$#source;$k++)
{
  $s_base = basename($source[$k]);
  $s_base =~ s/\.gz$//;
  $s_base =~ s/\.mnc$//;
  
  # 1 refine masks to avoid fitting in areas where information is missing
#  my $fast_xfm=run_elastix($elastix_fast,$source[$k],$target[$k],"$tmpdir/${s_base}_fast/");
  
#  do_cmd('itk_resample','--like',$source_mask[$k],$target_mask[$k],'--transform',$fast_xfm,"$tmpdir/${s_base}_fast/target_mask.mnc",'--labels','--invert');
  do_cmd('minccalc','-express','A[0]>0.5&&A[1]>0.5?1:0',$source_mask[$k],$target_mask[$k],"$tmpdir/${s_base}_fast/new_source_mask.mnc");
#  do_cmd('itk_resample','--like',$target_mask[$k],"$tmpdir/${s_base}_fast/new_source_mask.mnc",'--transform',$fast_xfm,"$tmpdir/${s_base}_fast/new_target_mask.mnc",'--labels');
  
  # zero out background
  
  do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',$source[$k],"$tmpdir/${s_base}_fast/new_source_mask.mnc","$tmpdir/${s_base}_fast/new_source.mnc");
  do_cmd('minccalc','-express','A[1]>0.5?A[0]:0',$target[$k],"$tmpdir/${s_base}_fast/new_source_mask.mnc","$tmpdir/${s_base}_fast/new_target.mnc");

  # 2 run real fit
  $tmp_xfm=run_elastix($elastix_par,"$tmpdir/${s_base}_fast/new_source.mnc","$tmpdir/${s_base}_fast/new_target.mnc","$tmpdir/$s_base/");#,"$tmpdir/${s_base}_fast/new_source_mask.mnc","$tmpdir/${s_base}_fast/new_target_mask.mnc"

  do_cmd('itk_resample','--like',$source_mask[$k],$target_mask[$k],'--transform',$tmp_xfm,"$tmpdir/$s_base/target_mask.mnc",'--labels','--invert');
  do_cmd('minccalc','-express','A[0]>0.5&&A[1]>0.5?1:0',$source_mask[$k],"$tmpdir/$s_base/target_mask.mnc","$tmpdir/$s_base/estimate_mask.mnc");

  if($opt{outdir})
  {
    do_cmd('cp',"$tmpdir/$s_base/TransformParameters.0.txt","$opt{outdir}/${s_base}_elx.txt");

    do_cmd('itk_resample',
          '--labels',
          '--like',$source_mask[$k],
          '--transform',$tmp_xfm,
          '--invert_transform',$bricks[$k],
          "$opt{outdir}/${s_base}_bricks.mnc");

    do_cmd('cp',"$tmpdir/$s_base/estimate_mask.mnc","$opt{outdir}/${s_base}_estimate_mask.mnc");
  }

  push @masks,"$tmpdir/$s_base/estimate_mask.mnc";
  push @grids,$tmp_grid;
}

regularize_grids(\@grids,\@masks,$opt{order},"$tmpdir/regularize.xfm",1);
cleanup_grids(\@grids) unless $opt{debug};

$prev_xfm = "$tmpdir/regularize.xfm";
$prev_grid = "$tmpdir/regularize_grid_0.mnc";

do_cmd('cp',"$tmpdir/regularize.par",$opt{par}) if $opt{par};
do_cmd('xfminvert',$prev_xfm,$outxfm);

sub do_cmd { 
   print STDOUT "@_\n" if $opt{verbose};
   if(!$opt{fake}){
      system(@_) == 0 or die "DIED: @_\n";
      }
   }
       
sub cleanup_xfms {
  my @xfms=@{$_[0]};
  for(my $i=0;$i<=$#xfms;$i++)
  {
    my $grid=$xfms[$i];
    $grid=~s/.xfm$/_grid_0.mnc/;
    do_cmd('rm','-f',$grid,$xfms[$i]);
  }
}

sub cleanup_grids {
  my @grids=@{$_[0]};
  for(my $i=0;$i<=$#grids;$i++)
  {
    do_cmd('rm','-f',$grids[$i]);
  }
}


sub regularize_grids {
  my ($_grid,$_mask,$order,$out,$invert,$out_par)=@_;
  my @grids=@{$_grid};
  my @masks=@{$_mask};
  die "XFMs and MASKs don't match!" if $#grids!=$#masks;
  
  my @args;
  if($opt{pca})
  {
    @args=('fit_harmonics_grids_regularize','--order',$opt{order},'--skip',2); # ANTS outputs grid @ 1mm
    push @args,'--cylindrical' if $opt{cyl};
    push @args,'--pca',$opt{pca};
    push @args,'--pcs',$opt{pcs} if $opt{pcs};
  } else {
    @args=($opt{cyl}?'c_fit_harmonics_grids':'fit_harmonics_grids','--order',$opt{order},'--skip',2);
    push(@args,"--limit") if $opt{limit};
    push(@args,"--keep",$opt{keep}) if $opt{keep};
    push(@args,'--iter',1) if $opt{keep}>0.99;
  }
  
  for(my $i=0;$i<=$#grids;$i++)
  {
    die "Can't find  grid file: $grids[$i]" unless -e $grids[$i];
    push (@args,$grids[$i],$masks[$i]);
  }
  
  push(@args,"$tmpdir/regularize.par",'--clobber');
  
  do_cmd(@args);
  @args=('par2xfm.pl',"$tmpdir/regularize.par",$out,'--clobber','--max',30,'--extent',400,'--step',4);
  push @args,'--noinvert'  if $invert;
  push @args,'--cylindric' if $opt{cyl};
  do_cmd(@args);
  do_cmd('cp',"$tmpdir/regularize.par",$out_par) if $out_par;
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}



sub run_elastix {
  my ($parameters,$source,$target,$out_dir,$source_mask,$target_mask)=@_;

  do_cmd('mkdir', '-p', $out_dir);

  open EL, ">$out_dir/elastix_parameters.txt" or die;
  print EL $parameters;
  close EL;

  my @args = ('elastix',
    '-f',  $source,
    '-m',  $target,
    '-out',  $out_dir, 
    '-p',  "$out_dir/elastix_parameters.txt",
    '-threads',4 );
 
  push (@args,'-fMask', $source_mask) if defined( $source_mask);
  push (@args,'-mMask',  $target_mask) if defined( $target_mask);

  do_cmd(@args);
  
  do_cmd('transformix', '-tp',  "$out_dir/TransformParameters.0.txt",
    '-def',  'all', '-out',  "$out_dir/");

  $tmp_grid = "$out_dir/deformationField.mnc";
  $tmp_xfm  = "$out_dir/deformationField.xfm";
  #   
  open XFM, ">$tmp_xfm" or die;
  print XFM "MNI Transform File\n";
  print XFM "Transform_Type = Linear;\nLinear_Transform =\n 1 0 0 0\n 0 1 0 0\n 0 0 1 0;\n";
  print XFM "Transform_Type = Grid_Transform;\n";
  print XFM "Displacement_Volume = deformationField.mnc;";
  close XFM;
  return $tmp_xfm;
}