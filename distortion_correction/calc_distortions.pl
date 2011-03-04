#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $coeff;
my $me=basename($0);
my $convert='mincconvert.static';
my $spacing=2;
my $output_tags;
my $output_model_tags;
my $order=3;
my $skip_grid=0;
my $align_xfm;

my $history=localtime() .">>> ".$me." ".join(' ',@ARGV);

GetOptions (    
	        "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "coeff=s"   => \$coeff,
          "spacing=f" => \$spacing,
          "tags=s"    => \$output_tags,
          "model_tags=s" => \$output_model_tags,
          "skip_grid"  => \$skip_grid,
          "align_xfm=s" =>\$align_xfm); 

die "Programm usage: $me <acr_phantom>.mnc <model_base_name> <output.xfm> [--clobber] [--verbose] [--coeff <coeffiecient_file> [--spacing <n>] [--tags <output_tags] [--align_xfm <xfm>] [--model_tags <model_tags>] [--skip_grid] \n" if $#ARGV<2;

my ($input,$model,$output)=@ARGV;
my $model_tags=$model;
$model_tags=~s/\.mnc|\.mnc\.gz/.tag/;

die "${output} exists!\n " if -e $output && !$clobber;
die "${model} doesn't exists!\n" if !-e $model;
die "${model_tags} doesn't exists!\n" if !-e $model_tags;
die "${output_tags} exists" if $output_tags && -e $output_tags && !$clobber;
die "${align_xfm} exists" if $align_xfm && -e $align_xfm && !$clobber;
die "${output_model_tags} exists" if $output_model_tags && -e $output_model_tags && !$clobber; 
my $output_grid=$output;
$output_grid=~s/\.xfm/_grid_0.mnc/;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

my $tmp_aligned="${tmpdir}/aligned.mnc";
my $tmp_aligned2="${tmpdir}/aligned2.mnc";
my $tmp_align_xfm="${tmpdir}/align.xfm";
my $tmp_align_xfm_i="${tmpdir}/align_i.xfm";
my $tmp_aligned_tags="${tmpdir}/aligned.tag";
my $tmp_native_tags="${tmpdir}/native.tag";
my $tmp_model_tags="${tmpdir}/model.tag";
my $tmp_model_native_tags="${tmpdir}/model_native.tag";
my $tmp_grid="${tmpdir}/grid.mnc";

#1. align input volume
do_cmd('acr_align.pl',$input,$model,$tmp_align_xfm);
do_cmd('mincresample',$input,'-transformation',$tmp_align_xfm,$tmp_aligned,'-like',$model);
do_cmd($convert,$tmp_aligned,$tmp_aligned2,'-2');

#2. measure volume
do_cmd('acr_measure',$tmp_aligned2,$tmp_aligned_tags);
#3 transform tags back into native space
do_cmd('xfminvert',$tmp_align_xfm,$tmp_align_xfm_i);
do_cmd('transformtags','-vol1','-transformation',$tmp_align_xfm_i,$model_tags,$tmp_model_native_tags);
do_cmd('transformtags','-vol1','-transformation',$tmp_align_xfm_i,$tmp_aligned_tags,$tmp_native_tags);
#4 calculate distortion
my @args=('calculate_distortions', $tmp_native_tags, $tmp_model_native_tags, $tmp_grid, '--order',$order,'--spacing',$spacing);
push @args,('--coeff',$coeff) if $coeff;
push @args,('--skip_grid') if $skip_grid;
do_cmd(@args);

#5 make xfm file
open  OF,">$output" or die "Can't open ${output} for writing!\n";
print OF "MNI Transform File\nTransform_Type = Grid_Transform;\n";
print OF "Invert_Flag = True;\n";
print OF "Displacement_Volume =";
print OF basename($output_grid);
print OF ";";
close OF;

#6 convert grid file
do_cmd($convert,$tmp_grid,$output_grid,'-clobber');

#7 modify history
do_cmd('minc_modify_header','-sinsert',':history='.$history,$output_grid);
#copy tags if needed
do_cmd('cp',$tmp_native_tags,$output_tags) if $output_tags;
do_cmd('cp',$tmp_model_native_tags,$output_model_tags) if $output_model_tags;
do_cmd('cp',$tmp_align_xfm,$align_xfm)     if $align_xfm;


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
        system(@_) == 0 or die "DIED: @_\n";
    }
}
