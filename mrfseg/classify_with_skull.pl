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
          "mask=s"    => \$mask,
          "separate=s" => \$separate
          );

die("Usage : $me <input_t1> <output_base> --mask mask [--clobber --verbose --separate <base>]\n") if $#ARGV<1;
my ($in,$out)=@ARGV;

die "Please specify mask (--mask <mask>)\n" unless $mask;
check_file($out) unless $clobber;


my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

#1 expand mask to include skull
do_cmd('itk_morph','--exp','D[4]',$mask,"$tmpdir/mask_d4.mnc");

#1 perform mrfsegmentation with mask 
do_cmd('gamixture',$in,"$tmpdir/mask_d4.mnc","$base/image_lowres_spec_prob_skull.txt","$tmpdir/fmm_skull",'-restarts',10);
do_cmd('mrfseg',$in,"$tmpdir/mask_d4.mnc","$base/image_lowres_spec_prob_skull.txt","$tmpdir/fmm_skull","$tmpdir/cls_skull.mnc");

#quick fix for data reordering in mrfseg
do_cmd('mincresample','-like',$mask,"$tmpdir/cls_skull.mnc","$tmpdir/cls_skull_.mnc",'-nearest');
do_cmd('mv',"$tmpdir/cls_skull_.mnc","$tmpdir/cls_skull.mnc");
#2 take skull&csf 
do_cmd('minccalc','-express','abs(A[0]-1)<0.5||abs(A[0]-4)<0.5',"$tmpdir/cls_skull.mnc","$tmpdir/mask_skull.mnc");
do_cmd('itk_morph','--exp','D[2] E[2]',"$tmpdir/mask_skull.mnc","$tmpdir/mask_skull_de.mnc");
do_cmd('mincmorph','-succ','GB[1:1:1:0]',"$tmpdir/mask_skull_de.mnc","$tmpdir/mask_skull_select.mnc");

#3 remove skull from the brain mask
do_cmd('mincresample','-like',$mask,"$tmpdir/mask_skull_select.mnc","$tmpdir/mask_skull_select_.mnc",'-nearest');
do_cmd('mv',"$tmpdir/mask_skull_select_.mnc","$tmpdir/mask_skull_select.mnc");
do_cmd('minccalc','-express','A[0]>0.5&&A[1]<0.5?1:0','-byte',$mask,"$tmpdir/mask_skull_select.mnc","$tmpdir/mask.mnc");

#4 second run of mrf segmentation
do_cmd('gamixture',$in,"$tmpdir/mask.mnc","$base/image_lowres_spec_prob.txt","$tmpdir/fmm",'-restarts',10);
do_cmd('mrfseg',$in,"$tmpdir/mask.mnc","$base/image_lowres_spec_prob.txt","$tmpdir/fmm","$tmpdir/cls_hard.mnc","$tmpdir/cls_pve.mnc");
#quick fix
do_cmd('mincresample','-like',$mask,"$tmpdir/cls_hard.mnc","$tmpdir/cls_hard_.mnc",'-nearest');
do_cmd('mincresample','-like',$mask,"$tmpdir/cls_pve.mnc","$tmpdir/cls_pve_.mnc",'-nearest');
do_cmd('mv',"$tmpdir/cls_hard_.mnc","$tmpdir/cls_hard.mnc");
do_cmd('mv',"$tmpdir/cls_pve_.mnc","$tmpdir/cls_pve.mnc");

#5 rename labels 
do_cmd('minccalc','-byte','-express','if(abs(A[0]-6)<0.5) 3 else if(abs(A[1]-3)<0.5) 4 else if(abs(A[1]-2)<0.5) 2 else 0',"$tmpdir/cls_pve.mnc","$tmpdir/cls_hard.mnc",$out,'-clobber');


if($separate)
{
  do_cmd('minccalc','-express','abs(A[0]-2)<0.5?1:0','-byte',$out,"${separate}_gm.mnc");
  do_cmd('minccalc','-express','abs(A[0]-4)<0.5?1:0','-byte',$out,"${separate}_wm.mnc");
  do_cmd('minccalc','-express','abs(A[0]-3)<0.5?1:0','-byte',$out,"${separate}_gmwm.mnc");
  do_cmd('minccalc','-express','A[0]>0.5?1:0',"$tmpdir/cls_hard.mnc","${separate}_brain.mnc");
}

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
