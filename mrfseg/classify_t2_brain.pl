#! /usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;


my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $mask;
my $base=dirname($0);
my $separate;

GetOptions (    
          "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          );

die("Usage : $me <input_t1> <output_pve> <output_mask>  [--clobber --verbose]\n") if $#ARGV<2;
my ($in,$out_pve,$out_msk)=@ARGV;

check_file($out_pve) unless $clobber;
check_file($out_msk) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

#1 perform mrfsegmentation with mask 
do_cmd('gamixture',$in,"default","$base/image_lowres_spec_prob_t2_2class.txt","$tmpdir/fmm_skull",'-restarts',10,"-parzensigma",2);
do_cmd('mrfseg',$in,"default","$base/image_lowres_spec_prob_t2_2class.txt","$tmpdir/fmm_skull","$tmpdir/cls_skull.mnc","$tmpdir/pve_skull.mnc");

#quick fix for data reordering in mrfseg
do_cmd('minccalc','-byte','-express','abs(A[0]-1)<0.5?1:0',"$tmpdir/pve_skull.mnc",$out_msk,'-clobber');
do_cmd('cp',"$tmpdir/pve_skull.mnc",$out_pve);


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
