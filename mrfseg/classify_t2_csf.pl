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
my $mask;
my $skull;

GetOptions (    
          "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "mask=s"    => \$mask,
          "skull=s"   => \$skull
          );

die("Usage : $me <input_t1> <output_csf> <output_pve> --mask <brain mask>  [--clobber --verbose]\n") if $#ARGV<1;
my ($in,$out,$out_pve)=@ARGV;

die "Provide --mask \n" unless $mask ;

check_file($out_pve) unless $clobber;
check_file($out) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
if($skull)
{
  do_cmd('minccalc','-express','A[0]>0.5&&A[1]<0.5?1:0',$mask,$skull,"$tmpdir/mask.mnc");
} else {
  do_cmd('minccalc','-express','A[0]>0.5?1:0',$mask,"$tmpdir/mask.mnc");
}
# todo: fill internal holes?
do_cmd('mincmorph','-succ','GB[1:1:1:0]',"$tmpdir/mask.mnc","$tmpdir/mask_s.mnc");
do_cmd('mv',"$tmpdir/mask_s.mnc","$tmpdir/mask.mnc");
#do_cmd('minccalc','-express','A[0]>0.5?0:1',"$tmpdir/mask_s.mnc","$tmpdir/mask_i.mnc");
#do_cmd('mincmorph','-succ','GB[1:1:0:1]',"$tmpdir/mask_i.mnc","$tmpdir/mask.mnc",'-clobber');

#fix a glitch in mincmorph
#my ($x,$y,$z)=split(/\n/,`mincinfo -dimlength xspace -dimlength yspace -dimlength zspace `);

#1 get a boundary layer
do_cmd('itk_morph','--exp','E[1]',"$tmpdir/mask.mnc","$tmpdir/e1.mnc");
do_cmd('itk_morph','--exp','E[4]',"$tmpdir/e1.mnc","$tmpdir/e5.mnc");

do_cmd('minccalc','-express','A[0]>0.5&&A[1]<0.5',"$tmpdir/mask.mnc","$tmpdir/e1.mnc","$tmpdir/b1.mnc");
do_cmd('minccalc','-express','A[0]>0.5&&A[1]<0.5',"$tmpdir/mask.mnc","$tmpdir/e5.mnc","$tmpdir/b5.mnc");
do_cmd('mincmorph','-succ','GB[1:1:1:0]',"$tmpdir/b5.mnc","$tmpdir/b5_.mnc");
do_cmd('mincmorph','-succ','GB[1:1:1:0]',"$tmpdir/b1.mnc","$tmpdir/b1_.mnc");
do_cmd('mv',"$tmpdir/b5_.mnc","$tmpdir/b5.mnc");
do_cmd('mv',"$tmpdir/b1_.mnc","$tmpdir/b1.mnc");

#1 perform mrfsegmentation with mask 
do_cmd('gamixture',$in,"$tmpdir/b5.mnc","$base/image_lowres_spec_prob_t2_2class.txt","$tmpdir/fmm_csf",'-restarts',10,"-parzensigma",2);
do_cmd('mrfseg',$in,"$tmpdir/b5.mnc","$base/image_lowres_spec_prob_t2_2class.txt","$tmpdir/fmm_csf","$tmpdir/cls_csf.mnc");

#quick fix for data reordering in mrfseg
do_cmd('minccalc','-byte','-express','abs(A[0]-2)<0.5||A[1]>0.5?1:0',"$tmpdir/cls_csf.mnc","$tmpdir/b1.mnc",$out,'-clobber');

#estimate mean and variance
my ($mean,$var)=split(/\n/,`mincstats -mean -var -q $in -mask $out -mask_binvalue 1 `);

#calculate partial volume map
do_cmd('minccalc','-express',"A[1]>0.5?100*exp(-(A[0]-$mean)*(A[0]-$mean)/(2*$var)):0",$in,"$tmpdir/mask.mnc",$out_pve,'-clobber');

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
