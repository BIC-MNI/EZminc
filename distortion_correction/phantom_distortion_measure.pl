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
my ($fit_xfm,$ideal,$core);
my $work_dir;
my $measure;
my $mask;
my $ex_r;
my $limit_linear=0;
my $keep=1.0;
my $cylindric=0;
my $acr;
my $init;
my $step_iterations;
my $measure;
my $debug;
my $adni;
my $out_roi;
my $only_roi;
my $dilate_roi=0;
my $pca;
my $pcs;
my $use_dd=0;
my $keep_tmp=0;

#additional parameters
GetOptions( 
          "verbose"   =>       \$verbose,
          "debug"     =>       \$debug,
          "clobber"   =>       \$clobber,
          "order=n"   =>       \$order,
          "work_dir=s"=>       \$work_dir,
          "model=s"   =>       \$model,
          "min_step=f"=>       \$min_step,
          "no_core_extract" => \$no_core_extract,
          "measure=s" =>       \$measure,
          "mask=s"      =>     \$mask,
          "exclude=f" =>       \$ex_r,
          "limit_linear" =>    \$limit_linear,
          "keep=f"    =>       \$keep,
          "cylindric" =>       \$cylindric,
          "acr=s"     =>       \$acr,
          "init=s"    =>       \$init,
          "step_iterations=n"=>\$step_iterations,
          "measure=s" =>       \$measure,
          "adni"      =>       \$adni,
          "out-roi=s" =>       \$out_roi,
          "only-roi"  =>       \$only_roi,
          "dd"        =>       \$use_dd,
          "dilate-roi=n" =>    \$dilate_roi,
          "pca=s"     =>       \$pca,
          "pcs=n"     =>       \$pcs,
          "keep-tmp"  =>       \$keep_tmp,
          );

die <<END 
 Usage: $me <scan_1> [<scan_2>...<scan_n>] <output.par> <output.xfm> 
  --model <model>       - ideal model 
[ 
  --mask  <mask file>   - mask the unwonted areas of the model (i.e edges)
  --order <n>           - approximation order  (default $order)
  --work_dir <dir>      - use this directory to keep intermediate files (useful for debugging)
  --min_step <n>        - minimal step size for the distortion approximation  (default $min_step) 
  --no_core_extract     - don't preprocess images  
  --limit_linear        - don't use
  --keep <n>            - perform LTS approximation with this fraction, when ==1 perform LSQ approxiamtion (default $keep)
  --cylindric           - assume cylindric (around Z axis) simmetricity of distortion field  
  --acr <real_acr>      - in case of ACR phantom specify the 'real' acr scan 
  --init <init_xfm>     - specify initial distortion estimation 
  --step_iterations <n> - specify number of iterations at each level of detail (used for debuggin mainly)
  --measure <output>    - perform residual deformation measurement at each iteration (long) 
  --debug               - for debugging
  --adni                - treat the phantom as adni 
  --out-roi <output>    - create ROI describing volume where distortions may be applied
  --only-roi            - a HACK to skip actual distortion correction field calculations 
  --pca <rotation.csv>  - Principal Components rotation matrix
  --pcs <n>             - Number of Principal components to use
  --dd                  - Use Diffeomorphic Demons registration instead of minctracc 
] 
END
  if $#ARGV<2; #number of arguments -1 

die "You need to specify a model (--model)!\n" unless $model && -e $model;

my ($scan,$output_par,$output_xfm);
$output_xfm=pop @ARGV;
$output_par=pop @ARGV;

my @scans=@ARGV;

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
}else{
  do_cmd('mkdir','-p',$work_dir);
}

if($ex_r) #making a spherical mask with specified radius
{
  $mask="$tmpdir/mask_sph.mnc";
  unless( -e  $mask)
  {
    do_cmd('make_phantom','-byte','-no_partial','-ellipse',
     '-width',$ex_r,$ex_r,$ex_r,
     '-center',0,0,0,
     '-start',-400,-400,-400,
     '-nelements',200,200,200,
     '-step',4,4,4,
     "$tmpdir/sphere.mnc");
    do_cmd('minccalc','-expression','1-A[0]',"$tmpdir/sphere.mnc",$mask);
  }
}

if(!$acr && !$adni)
{

  unless( -e "$tmpdir/skeleton.mnc") #creating a mask for ROI where fitting will happen
  {
    do_cmd('minccalc','-expres','A[0]>1?1:0',$model,"$tmpdir/skeleton.mnc") ;

    if($mask) {
      do_cmd('mincresample','-nearest','-like',"$tmpdir/skeleton.mnc",$mask,"$tmpdir/mask.mnc",'-q');
      do_cmd('minccalc','-expression','A[0]>0&&A[1]>0?1:0',"$tmpdir/skeleton.mnc","$tmpdir/mask.mnc","$tmpdir/skeleton_.mnc");
      do_cmd('mv',"$tmpdir/skeleton_.mnc","$tmpdir/skeleton.mnc");
    }
  }

  unless(-e "$tmpdir/fit.mnc")
  {
    if($mask)
    {
      #do_cmd('cp',$mask,"$tmpdir/fit.mnc");
      do_cmd('itk_morph','--threshold',1,'--exp','D[2] D[2] D[2] E[2]',$model,"$tmpdir/fit_1.mnc");
      do_cmd('mincresample','-nearest','-like',"$tmpdir/fit_1.mnc",$mask,"$tmpdir/mask.mnc",'-q','-clob');
      do_cmd('minccalc','-express','A[0]>0.5?A[1]:0',"$tmpdir/mask.mnc","$tmpdir/fit_1.mnc","$tmpdir/fit.mnc");
      do_cmd('rm','-f',"$tmpdir/fit_1.mnc");
    } else {
      do_cmd('itk_morph','--threshold',1,'--exp','D[2] D[2] D[2] E[2]',$model,"$tmpdir/fit.mnc");
    }
  }
} elsif( $adni) {
  if($mask) {
    do_cmd('cp',$mask,"$tmpdir/fit.mnc");
    do_cmd('cp',$mask,"$tmpdir/skeleton.mnc");
  } else {
    do_cmd('itk_morph','--threshold',1,'--exp','D[3]',$model,"$tmpdir/skeleton.mnc") 
        unless -e "$tmpdir/skeleton.mnc";
    do_cmd('cp',"$tmpdir/skeleton.mnc","$tmpdir/fit.mnc")
        unless -e "$tmpdir/fit.mnc";
  }
}

#now preparing files for each acquisition
my @args;
foreach $scan(@scans) {
  # extract core
  my $name=basename($scan,'.gz');
  unless(-e "$work_dir/core_${name}")
  {
    if($no_core_extract) # core was already extracted
    {
#       if($scan =~ /\.gz$/)
#       { 
#         do_cmd("gunzip -c $scan > $work_dir/core_${name}");
#       } else { 
        do_cmd('minccalc','-express','A[0]>20?A[0]:0',$scan,"$work_dir/core_${name}");
#       }
    } else {
      if($acr)
      {
        do_cmd('nu_correct',"-iter", 100, "-stop", 0.0001, "-fwhm", 0.1,$scan,"$tmpdir/nuc_${name}") unless -e "$tmpdir/nuc_${name}" ;
        my $t=`mincstats -biModalT -q $tmpdir/nuc_${name}`;
        $t=$t*1.0;
        my $max=`mincstats -q -pctT 98 $tmpdir/nuc_${name} -mask $tmpdir/nuc_${name} -mask_floor $t`;
        $max=$max*1.0;
        do_cmd('minccalc','-expression',"clamp(A[0]*100/$max,0,100)",
               "$tmpdir/nuc_${name}","$work_dir/core_${name}");
      } elsif( $adni ) {
        
        do_cmd('nu_correct',$scan,"-iter", 100, "-stop", 0.0001, "-fwhm", 0.1,"$tmpdir/n3_${name}") unless -e "$tmpdir/n3_${name}";
        do_cmd('make_phantom','-ellipse','-center',0,0,0,'-no_partial','-width',40,40,40,
                '-start','-100','-100','-100','-nelements',100,100,100,"$tmpdir/adni_center.mnc") 
            unless -e "$tmpdir/adni_center.mnc";
            
        my @com=split(/\s/,`mincstats -com -q -world_only $tmpdir/n3_${name}`);
        die "Can'f find COM of ADNI phantom \n" if $#com<2;
        do_cmd('param2xfm','-translation',$com[0],$com[1],$com[2],"$tmpdir/mv_${name}.xfm",'-clob');
        do_cmd('mincresample','-nearest',"$tmpdir/adni_center.mnc",
               '-like',"$tmpdir/n3_${name}",
               '-transform',"$tmpdir/mv_${name}.xfm",
               "$tmpdir/center_${name}.mnc",'-clob');
        
        my $threshold=`mincstats -q -mean $tmpdir/n3_${name} -mask $tmpdir/center_${name}.mnc -mask_binvalue 1`;
        chomp($threshold);
        $threshold*=0.7;
        do_cmd('minccalc','-express',"clamp(A[0]*100/$threshold,0,100)",
               "$tmpdir/n3_${name}","$work_dir/core_${name}") 
          unless -e "$work_dir/core_${name}";
          
      } else {
        # launch external programm for extracting 'core'
        do_cmd('lego_core_extract.pl',$scan,"$tmpdir/core_${name}",'--grayscale','--denoise','--nuc');
        do_cmd('minccalc','-express','A[0]>20?A[0]:0',"$tmpdir/core_${name}","$work_dir/core_${name}");
        do_cmd('rm','-f',"$tmpdir/core_${name}");
      }
    }
  }
  # align ideal representation to the scan
  if($acr && ! -e "$work_dir/align_${name}.xfm" ) #ADNI scan
  {
    # use predefined distortion correction to improve initial alignment
    if($init) {
      do_cmd('uniformize_minc.pl',"$work_dir/core_${name}","$tmpdir/init_${name}",
             '--transform',$init,'--resample','trilinear','--clobber') unless -e "$tmpdir/init_${name}";
      do_cmd('acr_align.pl',"$tmpdir/init_${name}",$acr,"$tmpdir/pre_align_${name}.xfm") unless -e "$tmpdir/pre_align_${name}.xfm";
    } else {
      do_cmd('acr_align.pl',"$tmpdir/core_${name}",$acr,"$tmpdir/pre_align_${name}.xfm") unless -e "$tmpdir/pre_align_${name}.xfm";
    }
    do_cmd('xfminvert',"$tmpdir/pre_align_${name}.xfm","$work_dir/align_${name}.xfm");
  }else{ #LEGO phantom scan
    if($init)
    {
      do_cmd('uniformize_minc.pl',"$tmpdir/core_${name}","$tmpdir/init_${name}",
             '--transform',$init,'--resample','trilinear','--clobber') unless -e "$tmpdir/init_${name}";

      do_cmd('bestlinreg.pl',$model,"$tmpdir/init_${name}",
             "$work_dir/align_${name}.xfm",'-lsq6') unless -e "$work_dir/align_${name}.xfm";

    } else {
     print "Manual xfm:$work_dir/mv_${name}.xfm\n";
     if(-e "$work_dir/mv_${name}.xfm") # allow specifiying initial transformation
     {
       print "Using manual xfm\n";
       do_cmd('bestlinreg_s',
              $model,"$work_dir/core_${name}",
              "$work_dir/align_${name}.xfm",'-lsq6',
              '-init_xfm',"$work_dir/mv_${name}.xfm") 
				unless -e "$work_dir/align_${name}.xfm";
     } else {
        do_cmd('bestlinreg.pl',$model,"$work_dir/core_${name}",
               "$work_dir/align_${name}.xfm",'-lsq6') 
        unless -e "$work_dir/align_${name}.xfm";
     }
    }
  }
  # create ideal representation

  $ENV{MINC_COMPRESS}=$minc_compress if $minc_compress;

  unless( -e "$work_dir/ideal_${name}") 
  { 
    do_cmd('mincresample',$model,"$tmpdir/ideal_${name}",'-like',
         "$work_dir/core_${name}",
         '-transform',"$work_dir/align_${name}.xfm") ;

    do_cmd('minccalc','-express','A[0]>20?A[0]:0',"$tmpdir/ideal_${name}","$work_dir/ideal_${name}");
  }

  # create mask
  unless(-e "$work_dir/fit_${name}")
  {
    if($acr)
    {
      die "ACR phantom processing needs --mask !\n" unless $mask;
      do_cmd('mincresample',$mask,"$work_dir/fit_${name}",
             '-like',"$work_dir/core_${name}",
             '-transform',"$work_dir/align_${name}.xfm",'-nearest');
    } else {

      do_cmd('mincresample',"$tmpdir/fit.mnc",
             "$work_dir/fit_${name}",'-like',"$work_dir/core_${name}",
             '-transform',"$work_dir/align_${name}.xfm",'-nearest');

      do_cmd('mincresample',"$tmpdir/skeleton.mnc","$tmpdir/estimate_${name}",
             '-like',"$work_dir/core_${name}",
             '-transform',"$work_dir/align_${name}.xfm",'-trilinear','-float');

    }
  }

  unless(-e "$work_dir/estimate_${name}" || $acr)
  {
      do_cmd('minccalc','-byte','-express','A[0]>0.5?1:0',"$tmpdir/estimate_${name}",
              "$work_dir/estimate_${name}");
  }

  if($acr) {
    push @args,"$work_dir/ideal_${name}","$work_dir/core_${name}","$work_dir/fit_${name}","$work_dir/fit_${name}";
  } else {
    push @args,"$work_dir/ideal_${name}","$work_dir/core_${name}","$work_dir/fit_${name}","$work_dir/estimate_${name}";
  }
}

unless( $only_roi )
{
  # calculate parameters
  if($use_dd)
  {
    @args=('phantomfit_DD.pl',@args,'-order',$order,'-clobber');
  } else {
    @args=('phantomfit.pl',@args,'-order',$order,'-clobber');

    push @args,'-weight',1;
    push @args,'-stiffness',0.4;
  }
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
  
  do_cmd(@args);
}

#create ROI
if($out_roi)
{
  my $output_field=$output_xfm;
  $output_field=~s/.xfm$/_grid_0.mnc/;

  my @masks;
  foreach $scan(@scans) {
    my $name=basename($scan,'.gz');
    if($dilate_roi>0)
    {
      do_cmd('itk_morph','--exp',"D[${dilate_roi}]",$mask,"$tmpdir/mask.mnc",'--clobber');
      $mask="$tmpdir/mask.mnc";
    }
    do_cmd('mincresample',$mask,'-like',$output_field,
           '-transform',"$work_dir/align_${name}.xfm",
           '-nearest',"$work_dir/mask_${name}.mnc",'-clobber');

    push @masks,"$work_dir/mask_${name}.mnc";
  }

  $ENV{MINC_COMPRESS}=$minc_compress if $minc_compress;
  do_cmd('mincmath','-byte','-max',@masks,$out_roi,'-clobber');
}

#exit 0 if $only_roi;



sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}

