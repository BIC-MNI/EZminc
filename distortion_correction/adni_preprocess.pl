#! /usr/bin/env perl

use strict; #become stricter

use File::Basename;            # for function basename
use File::Temp qw/ tempdir /;  # for temporary directory
use Getopt::Long;              # for parameters

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename( $0 ) ;

#additional parameters
GetOptions( 
          "verbose"   =>       \$verbose,
          "clobber"   =>       \$clobber,
          );
          
die <<END 
 Usage: $me <in> <out> 
[ 
  --clobber
  --verbose
] 
END
  if $#ARGV<1; #number of arguments -1 


my ($in,$out)=@ARGV;

check_file($out) unless $clobber;

my $minc_compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $minc_compress;

my  $tmpdir= &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );


do_cmd('mincnlm','-w',2,$in,"$tmpdir/denoised.mnc",'-mt',1);
do_cmd('nu_correct',"$tmpdir/denoised.mnc","-iter", 100, "-stop", 0.0001, "-fwhm", 0.1,"$tmpdir/n3.mnc");


do_cmd('make_phantom','-ellipse','-center',0,0,0,'-no_partial','-width',40,40,40,
      '-start','-100','-100','-100','-nelements',100,100,100,"$tmpdir/adni_center.mnc") 
  unless -e "$tmpdir/adni_center.mnc";
  
my @com=split(/\s/,`mincstats -com -q -world_only $tmpdir/n3.mnc`);
die "Can'f find COM of ADNI phantom \n" if $#com<2;

do_cmd('param2xfm','-translation',$com[0],$com[1],$com[2],"$tmpdir/mv.xfm",'-clob');
do_cmd('mincresample','-nearest',"$tmpdir/adni_center.mnc",
     '-like',"$tmpdir/n3.mnc",
     '-transform',"$tmpdir/mv.xfm",
     "$tmpdir/center.mnc",'-clob');

my $threshold=`mincstats -q -mean $tmpdir/n3.mnc -mask $tmpdir/center.mnc -mask_binvalue 1`;
chomp($threshold);
$threshold*=0.7;

do_cmd('minccalc','-express',"clamp(A[0]*100/$threshold,0,100)", "$tmpdir/n3.mnc",$out,'-clob') ;


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}

