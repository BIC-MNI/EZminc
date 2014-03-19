#! /usr/bin/env perl

############################# MNI Header #####################################
#@NAME       :  phantom_distortion_measure.pl
#@DESCRIPTION:  script to calculate distortions based on the phantom scans
#@COPYRIGHT  :
#              Vladimir S. Fonov  February, 2009
#              Montreal Neurological Institute, McGill University.
#              Permission to use, copy, modify, and distribute this
#              software and its documentation for any purpose and without
#              fee is hereby granted, provided that the above copyright
#              notice appear in all copies.  The author and McGill University
#              make no representations about the suitability of this
#              software for any purpose.  It is provided "as is" without
#              express or implied warranty.
###############################################################################


use strict; #become stricter

use File::Basename;            # for function basename
use File::Temp qw/ tempdir /;  # for temporary directory
use Getopt::Long;              # for parameters

my $fake=0;
my $verbose=0;
my $clobber=0;
my $order=5;
my $model;
my $min_step=2;
my $no_core_extract=0;
my $me=basename( $0 ) ;
my $mydir=dirname( $0 );
my ($fit_xfm,$ideal,$core);
my $work_dir;
my $measure;
my $mask;
my $ex_r;
my $limit_linear=0;
my $keep=1.0;
my $cylindric=0;
my $init;
my $step_iterations;
my $measure;
my $debug;
my $out_roi;
my $only_roi;
my $dilate_roi=0;
my $pca;
my $pcs;
my $use_dd=0;
my $use_ants=0;
my $keep_tmp=0;
my $use_elastix=0;
my $bricks;

#additional parameters
GetOptions( 
          "verbose"   =>       \$verbose,
          "debug"     =>       \$debug,
          "clobber"   =>       \$clobber,
          "order=n"   =>       \$order,
          "work_dir=s"=>       \$work_dir,
          "model=s"   =>       \$model,
          "min_step=f"=>       \$min_step,
          "measure=s" =>       \$measure,
          "mask=s"      =>     \$mask,
          "exclude=f" =>       \$ex_r,
          "limit_linear" =>    \$limit_linear,
          "keep=f"    =>       \$keep,
          "cylindric" =>       \$cylindric,
          "init=s"    =>       \$init,
          "step_iterations=n"=>\$step_iterations,
          "measure=s" =>       \$measure,
          "out-roi=s" =>       \$out_roi,
          "only-roi"  =>       \$only_roi,
          "dd"        =>       \$use_dd,
          "dilate-roi=n" =>    \$dilate_roi,
          "pca=s"     =>       \$pca,
          "pcs=n"     =>       \$pcs,
          "keep-tmp"  =>       \$keep_tmp,
          "ants"      =>       \$use_ants,
          "elastix"   =>       \$use_elastix, 
          "bricks=s" =>        \$bricks,
          );

die <<END 
Usage: $me <scan_1> <mask_1> [<scan_2> <mask_2>...<scan_n> <mask_n>] <output.par> <output.xfm> <output_measure.csv>
  --model <model>       - ideal model 
[ 
  --mask  <mask file>   - mask the unwonted areas of the model (i.e edges)
  --order <n>           - approximation order  (default $order)
  --work_dir <dir>      - use this directory to keep intermediate files (useful for debugging)
  --min_step <n>        - minimal step size for the distortion approximation  (default $min_step) 
  --limit_linear        - don't use
  --keep <n>            - perform LTS approximation with this fraction, when ==1 perform LSQ approxiamtion (default $keep)
  --cylindric           - assume cylindric (around Z axis) simmetricity of distortion field  
  --init <init_xfm>     - specify initial distortion estimation 
  --step_iterations <n> - specify number of iterations at each level of detail (used for debuggin mainly)
  --measure <output>    - perform residual deformation measurement at each iteration (long) 
  --debug               - for debugging
  --out-roi <output>    - create ROI describing volume where distortions may be applied
  --only-roi            - a HACK to skip actual distortion correction field calculations 
  --pca <rotation.csv>  - Principal Components rotation matrix
  --pcs <n>             - Number of Principal components to use
  --dd                  - Use Diffeomorphic Demons registration instead of minctracc 
  --ants                - Use mincANTS
  --elastix             - Use Elastix
  --bricks <bricks.mnc> - provide bricks file
] 
END
  if $#ARGV<4; #number of arguments -1 

die "You need to specify a model (--model) and mask (--mask)!\n" unless $model && -e $model && $mask && -e $mask;

my ($scan,$output_par,$output_xfm,$output_csv);

$output_csv=pop @ARGV;
$output_xfm=pop @ARGV;
$output_par=pop @ARGV;

my $i;
my @scans;
my @scan_masks;

for ($i=0;$i<=$#ARGV;$i+=2)
{  
  push(@scans, $ARGV[$i]);
  push(@scan_masks, $ARGV[$i+1]);
}

#check if the output exists
check_file($output_par) unless $clobber || $only_roi;
check_file($output_xfm) unless $clobber || $only_roi;
check_file($out_roi) if !$clobber && $out_roi;


#makes a temporary directory
my $tmpdir=&tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => !$keep_tmp );
my $minc_compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $minc_compress;

unless($work_dir)
{
  $work_dir=$tmpdir
} else {
  do_cmd('mkdir','-p',$work_dir);
}


#now preparing files for each acquisition
my @args;
for ($i=0;$i<=$#scans;$i+=1) {
  # extract core
  my $scan=$scans[$i];
  my $scan_mask=$scan_masks[$i];
  my $name=basename($scan,'.gz');
  
  unless(-e "$work_dir/${name}")
  {
    do_cmd('cp',$scan,"$work_dir/${name}");
  }
  
  # align ideal representation to the scan
  if( $init )
  {
    do_cmd('uniformize_minc.pl', "$tmpdir/${name}", "$tmpdir/init_${name}",
            '--transform',$init,'--resample','trilinear','--clobber') unless -e "$tmpdir/init_${name}";

    do_cmd('bestlinreg_s2', $model, "$tmpdir/init_${name}",
            "$work_dir/align_${name}.xfm",'-lsq6') unless -e "$work_dir/align_${name}.xfm";

  } else {
    if(-e "$work_dir/mv_${name}.xfm") # allow specifiying initial transformation
    {
      print "Using manual xfm:$work_dir/mv_${name}.xfm\n";
      do_cmd('bestlinreg_s2','-close',
            $model,"$work_dir/${name}",
            '-source_mask', $mask, 
            '-target_mask', $scan_mask, 
            "$work_dir/align_${name}.xfm",'-lsq6',
            '-init_xfm',"$work_dir/mv_${name}.xfm") 
        unless -e "$work_dir/align_${name}.xfm";
    } else {
      do_cmd('bestlinreg_s2',$model,"$work_dir/${name}",
              '-source_mask', $mask, 
              '-target_mask', $scan_mask, 
              "$work_dir/align_${name}.xfm",'-lsq6')
        unless -e "$work_dir/align_${name}.xfm";
    }
  }
  # create ideal representation
  $ENV{MINC_COMPRESS}=$minc_compress if $minc_compress;

  unless( -e "$work_dir/ideal_${name}") 
  { 
    do_cmd('mincresample',$model,"$work_dir/ideal_${name}",
           '-like',"$work_dir/${name}",
           '-transform',"$work_dir/align_${name}.xfm") ;
  }

  unless(-e "$work_dir/ideal_mask_${name}" )
  {
    do_cmd('itk_morph', '--threshold', '50', '--exp', 'D[4] E[2]', 
           "$work_dir/ideal_${name}", "$tmpdir/ideal_mask_${name}");
    
    do_cmd('mincresample',$mask,"$tmpdir/ideal_mask2_${name}",
          '-like',"$tmpdir/ideal_mask_${name}",
          '-transform',"$work_dir/align_${name}.xfm") ;

    do_cmd('minccalc', '-express', 'A[0]>0.5&&A[1]>0.5?1:0', '-byte', 
      "$tmpdir/ideal_mask_${name}", "$tmpdir/ideal_mask2_${name}",  "$work_dir/ideal_mask_${name}");

  }
  
  unless( -e "$work_dir/bricks_${name}") #consider using tfm_input_sampling here
  { 
    do_cmd('itk_resample',
           $bricks,"$work_dir/bricks_${name}",
           '--like',"$work_dir/${name}",
           '--transform',"$work_dir/align_${name}.xfm",
           '--labels') ;
  }

  # create mask
  unless(-e "$work_dir/mask_${name}")
  {
      
      do_cmd('itk_morph', '--threshold', '50', '--exp', 'D[4] E[2]', "$work_dir/${name}", "$tmpdir/mask_${name}");
      do_cmd('mincresample','-nearest','-like', "$tmpdir/mask_${name}", $scan_mask, "$tmpdir/mask2_${name}");
      do_cmd('minccalc', '-express', 'A[0]>0.5&&A[1]>0.5?1:0', '-byte', 
        "$tmpdir/mask_${name}", "$tmpdir/mask2_${name}", "$work_dir/mask_${name}");
  }

  
  push @args,"$work_dir/ideal_${name}","$work_dir/${name}","$work_dir/ideal_mask_${name}","$work_dir/mask_${name}","$work_dir/bricks_${name}";
}

do_cmd('mkdir','-p',"$work_dir/work") if $debug;

unless( $only_roi )
{
  # calculate parameters
  @args=("$mydir/phantomfit_elastix.pl",@args,'-order',$order,'-clobber');

  push @args,'-work_dir',"$work_dir/work"  if $debug ;
  push @args,'-cylindric'        if $cylindric;
  push @args,'-keep',$keep       if $keep;
  push @args,"-measure",$measure if $measure;
  push @args,"-limit"            if $limit_linear;
  push @args,"-init",$init       if $init;
  push @args,'-par',$output_par,'-min_step',$min_step,$output_xfm;
  push @args,'-measure',$measure if $measure;
  push @args,'-step_iterations',$step_iterations if $step_iterations;
  push @args,'-debug' if $debug;
  push @args,'-pca',$pca if $pca;
  push @args,'-pcs',$pcs if $pcs;
  push @args,'-outdir',$work_dir;

  do_cmd(@args);
}

#create ROI
if($out_roi)
{
  my $output_field=$output_xfm;
  $output_field=~s/.xfm$/_grid_0.mnc/;

  my @masks;
  for ($i=0;$i<=$#scans;$i+=1) {
# extract core
    my $scan=$scans[$i];
    my $scan_mask=$scan_masks[$i];
    my $name=basename($scan,'.gz');

    do_cmd('mincresample',$scan_mask,'-like',$output_field,
           '-nearest',"$work_dir/mask_lr_${name}.mnc",'-clobber');

    push @masks,"$work_dir/mask_lr_${name}.mnc";
  }

  $ENV{MINC_COMPRESS}=$minc_compress if $minc_compress;
  do_cmd('mincmath','-byte','-max',@masks,$out_roi,'-clobber');
}

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}

