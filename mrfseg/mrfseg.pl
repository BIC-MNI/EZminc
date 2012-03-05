#! /usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempdir /;
use Getopt::Long;


my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $mask;
my $base=dirname($0);
my $atlas;
my $lowres;

GetOptions( 
          "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "mask=s"    => \$mask,
          "atlas=s"   => \$atlas,
          "lowres"    => \$lowres
          
          );

die("Usage : $me <input_t1> <output.mask> [--mask mask --clobber --verbose --atlas <atlas spec> --lowres]\n") if $#ARGV<1;
my ($in,$out)=@ARGV;


unless($atlas)
{
  unless($lowres)
  {
    $atlas="$base/icbm_image_spec_prob.txt" ;
  }else {
    $atlas="$base/image_lowres_spec_prob.txt";
  } 
}
$mask="default" unless $mask;

check_file($out) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

unless($lowres)
{
  do_cmd('gamixture',$in,$mask,$atlas,"$tmpdir/fmm",'-restarts',10);
} else {
  do_cmd('gamixture',$in,$mask,$atlas,"$tmpdir/fmm",'-restarts',10,'-parzensigma',2);
}
do_cmd('mrfseg',$in,$mask,$atlas,"$tmpdir/fmm",$out);

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
}
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
