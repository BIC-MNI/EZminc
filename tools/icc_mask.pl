#! /usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempdir /;
use Getopt::Long;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $dilate=0;
my $model;
my $icc_model;

GetOptions (    
          "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "model=s"   => \$model,
          "icc-model=s" => \$icc_model
          );
          
die "Usage: $me <scan_in> <brain_mask_out> --model <T1w model> --icc-model <icc mask model> \n[--verbose\n --clobber\n ]\n" if $#ARGV<1;

my ($in,$out)=@ARGV;

die "Specify model and icc-model \n" unless $model && $icc_model;

check_file($out) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

do_cmd('nlfit_s','-level',4,$in,$model,"$tmpdir/reg.xfm");
do_cmd('mincresample','-nearest','-like',$in,$icc_model,'-transform',"$tmpdir/reg.xfm",'-invert_transformation',$out,'-clob');


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
}
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
