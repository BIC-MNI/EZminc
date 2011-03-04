#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
#use POSIX qw(floor);

my $fake=0;
my $verbose=0;
my $clobber=0;
my $keep_tmp=0;
my $me=basename($0);
my $grayscale;
my $nuc;
my $denoise;
my $mri_3t;
my $nuc_iter=1;

GetOptions (    
        "verbose"   => \$verbose,
        "clobber"   => \$clobber,
        "grayscale" => \$grayscale,
        "nuc"       => \$nuc,
        "denoise"   => \$denoise,
        '3t'        => \$mri_3t,
        'nuc-iter=n'=> \$nuc_iter
          );

die "Usage: $me <input> <output> [--verbose --clobber --grayscale --nuc --denoise --3t --nuc-iter <n>] \n" if $#ARGV<1;

my ($input,$output)=@ARGV;

check_file($output) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => !$keep_tmp );
#add black padding



my $comp=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $comp;

#0. preprocess
do_cmd('mincreshape',$input,'+direction',
       '-dimsize','xspace=-1',
       '-dimsize','yspace=-1',
       '-dimsize','zspace=-1',
       '-dimorder','zspace,yspace,xspace',
       "$tmpdir/input.mnc",'-float');

my ($spx,$spy,$spz)=split(/\n/,`mincinfo -attvalue xspace:spacing -attvalue yspace:spacing -attvalue zspace:spacing $tmpdir/input.mnc`);

do_cmd('minc_modify_header','-sinsert','xspace:spacing=regular__',"$tmpdir/input.mnc") if $spx=='irregular';
do_cmd('minc_modify_header','-sinsert','yspace:spacing=regular__',"$tmpdir/input.mnc") if $spy=='irregular';
do_cmd('minc_modify_header','-sinsert','zspace:spacing=regular__',"$tmpdir/input.mnc") if $spz=='irregular';
$input="$tmpdir/input.mnc";

my ($xspace,$yspace,$zspace,$xstart,$ystart,$zstart,$xstep,$ystep,$zstep)= split(/\n/, `mincinfo  -dimlength xspace -dimlength yspace -dimlength zspace -attvalue xspace:start -attvalue yspace:start -attvalue zspace:start -attvalue xspace:step -attvalue yspace:step -attvalue zspace:step $input`);
$xspace+=20;
$yspace+=20;
$zspace+=20;

if($denoise)
{
  do_cmd('mincnlm',$input,"$tmpdir/denoised.mnc");
  $input="$tmpdir/denoised.mnc";
}

my @args;

unless($grayscale)
{
 nu_correct($input,$nuc_iter,"$tmpdir/corrected.mnc");
 $input="$tmpdir/corrected.mnc";
 
 do_cmd('mincreshape','-dimrange',"xspace=-10,$xspace",'-dimrange',"yspace=-10,$yspace",'-dimrange',"zspace=-10,$zspace",$input,"$tmpdir/padded.mnc");
 do_cmd('mv',"$tmpdir/padded.mnc","$tmpdir/corrected.mnc");
 do_cmd('itk_morph','--bimodal','--exp','D[4] E[5]',"$tmpdir/corrected.mnc","$tmpdir/mask.mnc");
 do_cmd('mincresample','-nearest',"$tmpdir/corrected.mnc",'-like',"$tmpdir/mask.mnc","$tmpdir/scan.mnc");
 my $threshold=`mincstats -mask $tmpdir/mask.mnc -mask_binvalue 1 -q -pctT 50 $tmpdir/scan.mnc`;
 chomp($threshold);
 do_cmd('minccalc','-clobber','-expression',"A[0]?clamp(100*($threshold-A[1])/($threshold),0,100):0",
        "$tmpdir/mask.mnc","$tmpdir/scan.mnc","$tmpdir/corrected.mnc");
} else {
  
 if($nuc) {
  nu_correct($input,$nuc_iter,"$tmpdir/nuc.mnc");
  $input="$tmpdir/nuc.mnc";
 }
 do_cmd('mincreshape','-dimrange',"xspace=-10,$xspace",'-dimrange',"yspace=-10,$yspace",'-dimrange',"zspace=-10,$zspace",$input,"$tmpdir/corrected.mnc",'-clobber');
 do_cmd('itk_g_morph','--exp','D[3]',"$tmpdir/corrected.mnc","$tmpdir/corrected_d3.mnc");
 do_cmd('itk_g_morph','--exp','E[3]',"$tmpdir/corrected_d3.mnc","$tmpdir/mask.mnc");
 do_cmd('mincresample','-nearest',"$tmpdir/corrected.mnc",'-like',"$tmpdir/mask.mnc","$tmpdir/scan.mnc");
 my $pct10=`mincstats -q -biModalT $tmpdir/input.mnc`; # -pctT 10
 chomp($pct10);
 do_cmd('minccalc','-expression',"A[2]>${pct10}&&(A[0]-A[1])>10?(A[0]-A[1])/A[2]:0","$tmpdir/mask.mnc","$tmpdir/scan.mnc","$tmpdir/corrected_d3.mnc","$tmpdir/diff.mnc");
 my $threshold=`mincstats -q -pctT 99 $tmpdir/diff.mnc`;
 chomp($threshold);
 $threshold/=100.0;
 do_cmd('minccalc','-clobber','-expression',"clamp(A[0]/${threshold},0,100)",
        "$tmpdir/diff.mnc","$tmpdir/corrected.mnc",'-byte');
}


$xspace-=20;
$yspace-=20;
$zspace-=20;

$ENV{MINC_COMPRESS}=$comp if $comp;
do_cmd('mincreshape','-dimrange',"xspace=10,$xspace",
                     '-dimrange',"yspace=10,$yspace",
                     '-dimrange',"zspace=10,$zspace",
                     "$tmpdir/corrected.mnc",$output,'-clobber');


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}


sub nu_correct {
 my ($input,$nuc_iter,$output)=@_;

 my $iter=0;
 while($iter<$nuc_iter)
 {
  @args=("nu_correct", "-clobber", "-iter", 100, "-stop", 0.00001, "-fwhm", 0.1,$input,  "$tmpdir/corrected_$iter.mnc",'-clobber');
  push @args,'-distance',50 if $mri_3t;
  if($iter>0)
  {
    do_cmd('itk_morph','--bimodal','--exp','E[1]',$input,"$tmpdir/mask_$iter.mnc");
    push(@args,'-mask',"$tmpdir/mask_$iter.mnc");
  }
  do_cmd(@args);
  $input="$tmpdir/corrected_$iter.mnc";
  $iter+=1;
 }
 do_cmd('mv',$input,$output);
}
