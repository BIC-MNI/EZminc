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
my $correct;

GetOptions (    
	        "verbose"    => \$verbose,
          "clobber"    => \$clobber,
          "mask=s"     => \$mask,
          "correct"    => \$correct,
          );

die("Usage : $me <input_t1> <output.mask> --mask mask [--clobber --correct  --verbose]\n") if $#ARGV<1;
my ($in,$out)=@ARGV;

die "Please specify mask (--mask <mask>)\n" unless $mask;
check_file($out) unless $clobber;


my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

#1 expand mask to include skull
do_cmd('itk_morph','--exp','D[4]',$mask,"$tmpdir/mask_d4.mnc");
do_cmd('itk_morph','--exp','E[4]',$mask,"$tmpdir/mask_e4.mnc");
do_cmd('minccalc','-express','A[0]>0.5&&A[1]<0.5?1:0','-byte',$mask,"$tmpdir/mask_e4.mnc","$tmpdir/mask_b8.mnc");

#1 perform mrfsegmentation with mask 
do_cmd('gamixture',$in,"$tmpdir/mask_d4.mnc","$base/image_lowres_spec_prob_skull.txt","$tmpdir/fmm_skull",'-restarts',10,"-parzensigma",2);
do_cmd('mrfseg',$in,"$tmpdir/mask_d4.mnc","$base/image_lowres_spec_prob_skull.txt","$tmpdir/fmm_skull","$tmpdir/cls_skull.mnc");

#quick fix for data reordering in mrfseg
do_cmd('mincresample','-like',$mask,"$tmpdir/cls_skull.mnc","$tmpdir/cls_skull_.mnc",'-nearest');
do_cmd('mv',"$tmpdir/cls_skull_.mnc","$tmpdir/cls_skull.mnc");
#2 take skull&csf 
do_cmd('minccalc','-express','abs(A[0]-4)<0.5||abs(A[0]-1)<0.5?1:0','-byte',"$tmpdir/cls_skull.mnc","$tmpdir/mask_skull.mnc");#||abs(A[0]-1)<0.5 for csf
do_cmd('itk_morph','--exp','D[3] E[3]',"$tmpdir/mask_skull.mnc","$tmpdir/mask_skull_de.mnc");
do_cmd('mincmorph','-succ','GB[1:1:1:0]',"$tmpdir/mask_skull_de.mnc","$tmpdir/mask_skull_select.mnc");

#3 remove skull from the brain mask
do_cmd('mincresample','-like',$mask,"$tmpdir/mask_skull_select.mnc","$tmpdir/skull_.mnc",'-nearest','-clobber');
if($correct)
{
  do_cmd('minccalc','-express','A[0]>0.5&&A[1]>0.5?1:0','-byte',"$tmpdir/skull_.mnc","$tmpdir/mask_b8.mnc","$tmpdir/mask_skull.mnc",'-clobber');
  do_cmd('minccalc','-express','A[0]>0.5&&A[1]<0.5?1:0','-byte',$mask,"$tmpdir/mask_skull.mnc","$tmpdir/brain.mnc",'-clobber');
  do_cmd('itk_morph','--exp','E[3] D[3]',"$tmpdir/brain.mnc","$tmpdir/brain_.mnc");
  do_cmd('mincmorph','-succ','GB[1:1:1:0]',"$tmpdir/brain_.mnc",$out,'-clobber');
} else {
  do_cmd('minccalc','-express','A[0]>0.5&&A[1]>0.5?1:0','-byte',"$tmpdir/skull_.mnc","$tmpdir/mask_b8.mnc",$out,'-clobber');
}
#do_cmd('mv',"$tmpdir/mask_skull_select_.mnc",$out);
#do_cmd('minccalc','-express','A[0]>0.5&&A[1]>0.5?1:0','-byte',"$tmpdir/mask_b4.mnc","$tmpdir/mask_skull_select.mnc",$out,'-clobber');


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
