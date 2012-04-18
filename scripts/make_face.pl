#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $fwhm=2;
my $keep_tmp=0;
my $rotate;
my $level=1.0;
my $smooth;
my $stx_xfm;

GetOptions (    
	        "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "fwhm=f"    => \$fwhm,
          "rotate=f"    => \$rotate,
          "threshold=f" => \$level,
          "smooth"    => \$smooth,
          'stx=s'     => \$stx_xfm
          ); 

die "Usage: $me <scan_in> <face_out> [--verbose --clobber --fwhm <f> --rotate <n> --threshold <f> --stx <xfm>] \n" if $#ARGV<1;

my ($scan,$output)=@ARGV;

check_file($output) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => !$keep_tmp );

delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};
#do_cmd('nu_correct', $scan,"$tmpdir/corr.mnc");
#blur

if($rotate)
{
  do_cmd('param2xfm','-rotations',0,0,$rotate,"$tmpdir/rotate.xfm");
}

if($stx_xfm) #make unscaled
{
  my $scale=`xfm2param $stx_xfm|fgrep scale`;
  chomp($scale);
  print $scale,"\n";
  do_cmd("param2xfm $scale $tmpdir/scale.xfm");
  do_cmd('xfminvert',"$tmpdir/scale.xfm","$tmpdir/unscale.xfm");
  do_cmd('xfmconcat',$stx_xfm,"$tmpdir/unscale.xfm","$tmpdir/stx.xfm");
  if($rotate)
  {
    do_cmd('xfmconcat',"$tmpdir/stx.xfm","$tmpdir/rotate.xfm","$tmpdir/rstx.xfm");
    do_cmd('cp',"$tmpdir/rstx.xfm","$tmpdir/rotate.xfm");
  } else {
    do_cmd('cp',"$tmpdir/stx.xfm","$tmpdir/rotate.xfm");
  }
} 

if($rotate||$stx_xfm)
{
  do_cmd('mincresample','-tfm_input_sampling',$scan,"$tmpdir/scan_.mnc",'-transformation',"$tmpdir/rotate.xfm");
  $scan="$tmpdir/scan_.mnc";
}

do_cmd('uniformize_minc.pl',$scan,"$tmpdir/resample.mnc",'--resample','trilinear','--step',1,'--datatype','byte');

my $threshold=`mincstats -q -biModalT $tmpdir/resample.mnc`;
chomp($threshold);
$threshold*=$level;
do_cmd('mincmorph','-successive',"B[$threshold:1000000:1:0]GB[1:1:1:0]","$tmpdir/resample.mnc","$tmpdir/scan.mnc");
#do_cmd('minccalc','-expression',"A[0]>$threshold?1:0",'-byte',"$tmpdir/resample.mnc","$tmpdir/scan.mnc");
$threshold=0.2;
do_cmd('rm','-f',"$tmpdir/resample.mnc") unless $keep_tmp;

do_cmd('mincblur',"$tmpdir/scan.mnc",'-fwhm',$fwhm,"$tmpdir/scan");
do_cmd('rm','-f',"$tmpdir/scan.mnc") unless $keep_tmp;

#do_cmd('itk_morph',"$tmpdir/resample.mnc","$tmpdir/scan_blur.mnc",'--bimodal','--exp','M[4]');
#marching cubes
do_cmd('marching_cubes',"$tmpdir/scan_blur.mnc","$tmpdir/scan.obj",$threshold);

#triangulate
do_cmd('triangulate_polygons',"$tmpdir/scan.obj","$tmpdir/scan2.obj");
#smooth surface
if($smooth)
{
  do_cmd('sm',"$tmpdir/scan2.obj","$tmpdir/scan3.obj",1000000,10);
  do_cmd('mv',"$tmpdir/scan3.obj","$tmpdir/scan2.obj");
}
#raytrace
do_cmd('ray_trace',"$tmpdir/scan2.obj",'-output',"$tmpdir/scan2.rgb",'-size',500,500,'-crop','-front');
#convert
do_cmd('convert',"$tmpdir/scan2.rgb",$output);

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}
sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
