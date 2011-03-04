#! /usr/bin/env perl

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
          "adni"      =>       \$adni
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
] 
END
  if $#ARGV<2; #number of arguments -1 

die "You need to specify a model (--model)!\n" unless $model && -e $model;

my ($scan,$output_par,$output_xfm);
$output_xfm=pop @ARGV;
$output_par=pop @ARGV;

my @scans=@ARGV;

#check if the output exists
check_file($output_par) unless $clobber;
check_file($output_xfm) unless $clobber;

#makes a temporary directory
my $tmpdir;

my $minc_compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $minc_compress;


unless($work_dir)
{
  $tmpdir= &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
}else{
  $tmpdir=$work_dir;
  do_cmd('mkdir','-p',$work_dir);
}

if($ex_r)
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

if(!$acr && ! $adni)
{
  unless( -e "$tmpdir/skeleton.mnc")
  {
    do_cmd('itk_morph','--threshold',1,'--exp','D[1]',$model,"$tmpdir/skeleton.mnc") ;
    if($mask) {
      do_cmd('minccalc','-expression','A[0]>0&&A[1]>0?1:0',"$tmpdir/skeleton.mnc",$mask,"$tmpdir/skeleton_.mnc");
      do_cmd('mv',"$tmpdir/skeleton_.mnc","$tmpdir/skeleton.mnc");
    }
  }
  unless(-e "$tmpdir/fit.mnc")
  {
    do_cmd('itk_morph','--threshold',1,'--exp','D[3]',$model,"$tmpdir/fit.mnc");
    if($mask) {
      do_cmd('minccalc','-expression','A[0]>0&&A[1]>0?1:0',"$tmpdir/fit.mnc",$mask,"$tmpdir/fit_.mnc");
      do_cmd('mv',"$tmpdir/fit_.mnc","$tmpdir/fit.mnc");
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

my @args;
foreach $scan(@scans) {
  # extract core
  my $name=basename($scan,'.gz');
  unless(-e "$tmpdir/core_${name}")
  {
    if($no_core_extract)
    {
      if($scan =~ /\.gz$/)
      { 
        do_cmd("gunzip -c $scan > $tmpdir/core_${name}");
      } else {
        do_cmd('cp',$scan,"$tmpdir/core_${name}");
      }
    } else {
      if($acr)
      {
        do_cmd('nu_correct',"-iter", 100, "-stop", 0.0001, "-fwhm", 0.1,$scan,"$tmpdir/nuc_${name}") unless -e "$tmpdir/nuc_${name}" ;
        my $t=`mincstats -biModalT -q $tmpdir/nuc_${name}`;
        $t=$t*1.0;
        my $max=`mincstats -q -pctT 98 $tmpdir/nuc_${name} -mask $tmpdir/nuc_${name} -mask_floor $t`;
        $max=$max*1.0;
        do_cmd('minccalc','-expression',"clamp(A[0]*100/$max,0,100)",
               "$tmpdir/nuc_${name}","$tmpdir/core_${name}");
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
               "$tmpdir/n3_${name}","$tmpdir/core_${name}") 
          unless -e "$tmpdir/core_${name}";
          
      } else {
        do_cmd('lego_core_extract.pl',$scan,"$tmpdir/core_${name}",'--grayscale');
      }
    }
  }
  # align ideal
  if($acr && ! -e "$tmpdir/align_${name}.xfm" )
  {
    if($init) {
      do_cmd('uniformize_minc.pl',"$tmpdir/core_${name}","$tmpdir/init_${name}",'--transform',$init,'--resample','trilinear','--clobber') unless -e "$tmpdir/init_${name}";
      do_cmd('acr_align.pl',"$tmpdir/init_${name}",$acr,"$tmpdir/pre_align_${name}.xfm") unless -e "$tmpdir/pre_align_${name}.xfm";
    } else {
      do_cmd('acr_align.pl',"$tmpdir/core_${name}",$acr,"$tmpdir/pre_align_${name}.xfm") unless -e "$tmpdir/pre_align_${name}.xfm";
    }
    do_cmd('xfminvert',"$tmpdir/pre_align_${name}.xfm","$tmpdir/align_${name}.xfm");
  }else{
    if($init)
    {
      do_cmd('uniformize_minc.pl',"$tmpdir/core_${name}","$tmpdir/init_${name}",'--transform',$init,'--resample','trilinear','--clobber') unless -e "$tmpdir/init_${name}";
      do_cmd('bestlinreg.pl',$model,"$tmpdir/init_${name}","$tmpdir/align_${name}.xfm",'-lsq6') unless -e "$tmpdir/align_${name}.xfm";
    } else {
     if(-e "$tmpdir/mv_${name}.xfm")
     {
       do_cmd('bestlinreg.pl',
              $model,"$tmpdir/core_${name}",
              "$tmpdir/align_${name}.xfm",'-lsq6',
              '-init_xfm',"$tmpdir/mv_${name}.xfm",'-noresample') 
				unless -e "$tmpdir/align_${name}.xfm";
     } else {
        do_cmd('bestlinreg.pl',$model,"$tmpdir/core_${name}",
               "$tmpdir/align_${name}.xfm",'-lsq6') 
        unless -e "$tmpdir/align_${name}.xfm";
     }
    }
  }
  # create ideal representation
  do_cmd('mincresample',$model,"$tmpdir/ideal_${name}",'-like',
         "$tmpdir/core_${name}",
         '-transform',"$tmpdir/align_${name}.xfm") unless -e "$tmpdir/ideal_${name}";
  # create mask
  unless(-e "$tmpdir/fit_${name}")
  {
    if($acr)
    {
      die "ACR phantom processing needs --mask !\n" unless $mask;
      do_cmd('mincresample',$mask,"$tmpdir/fit_${name}",'-like',"$tmpdir/core_${name}",'-transform',"$tmpdir/align_${name}.xfm",'-nearest');
    } else {
      do_cmd('mincresample',"$tmpdir/fit.mnc","$tmpdir/fit_${name}",'-like',"$tmpdir/core_${name}",'-transform',"$tmpdir/align_${name}.xfm",'-nearest');
      do_cmd('mincresample',"$tmpdir/skeleton.mnc","$tmpdir/estimate_${name}",'-like',"$tmpdir/core_${name}",'-transform',"$tmpdir/align_${name}.xfm",'-nearest');
    }
  }
  if($acr) {
    push @args,"$tmpdir/ideal_${name}","$tmpdir/core_${name}","$tmpdir/fit_${name}","$tmpdir/fit_${name}";
  } else {
    push @args,"$tmpdir/ideal_${name}","$tmpdir/core_${name}","$tmpdir/fit_${name}","$tmpdir/estimate_${name}";
  }
}

# calculate parameters
@args=('phantomfit.pl',@args,'-order',$order,'-clobber');
push @args,'-cylindric'        if $cylindric;
push @args,'-keep',$keep       if $keep;
push @args,"-measure",$measure if $measure;
push @args,"-limit"            if $limit_linear;
push @args,"-init",$init       if $init;
push @args,'-par',$output_par,'-min_step',$min_step,$output_xfm;
push @args,'-measure',$measure if $measure;
push @args,'-step_iterations',$step_iterations if $step_iterations;
push @args,'-work_dir',"$work_dir/pf" if $work_dir;
push @args,'-debug' if $debug;

#test
push @args,'-weight',1;
push @args,'-stiffness',0;


do_cmd(@args);

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}

