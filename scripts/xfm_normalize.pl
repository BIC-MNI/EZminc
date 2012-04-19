#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempdir /;
use Getopt::Long;
use POSIX qw(ceil);

my $fake=0;
my $verbose=0;
my $clobber=0;
my $model;
my $step;
my $me=basename($0);
my $exact;
my $invert;

GetOptions (    
        "verbose"   => \$verbose,
        "clobber"   => \$clobber,
        "like=s"    => \$model,
        "step=f"    => \$step,
        "exact"     => \$exact,
        "invert"    => \$invert);

die "Usage: $me <xfm_input|grid_input> <xfm_output> [--verbose --clobber --like <minc> --step <n> --exact --invert] \n" if $#ARGV<1;
my ($in,$out)=@ARGV;

$model=$in if $in =~ /.mnc$|.mnc.gz$/ && !$model;

die "You have to use --like <sample> or provide input grid file \n" unless $model;


my $out_grid=$out;
$out_grid=~s/\.xfm$//;
$out_grid.='_grid_0.mnc';

check_file($out) unless $clobber;
check_file($out_grid) unless $clobber;

my $tmpdir;


my ($xspace,$yspace,$zspace,$xstart,$ystart,$zstart,$xstep,$ystep,$zstep)= split(/\n/, `mincinfo  -dimlength xspace -dimlength yspace -dimlength zspace -attvalue xspace:start -attvalue yspace:start -attvalue zspace:start -attvalue xspace:step -attvalue yspace:step -attvalue zspace:step $model`);

$xspace=$xspace*1;
$yspace=$yspace*1;
$zspace=$zspace*1;

#I hate perl!
$xstart=$xstart*1.0;
$ystart=$ystart*1.0;
$zstart=$zstart*1.0;
#I hate perl!
$xstep=$xstep*1.0;
$ystep=$ystep*1.0;
$zstep=$zstep*1.0;

if($xstep<0.0)
{
  $xstart+=$xstep*($xspace-1);
  $xstep= -$xstep;
}

if($ystep<0.0)
{
  $ystart+=$ystep*($yspace-1);
  $ystep= -$ystep;
}

if($zstep<0.0)
{
  $zstart+=$zstep*($zspace-1);
  $zstep= -$zstep;
}

my $xlen=$xstep*$xspace;
my $ylen=$ystep*$yspace;
my $zlen=$zstep*$zspace;

if($step)
{
  $xstep=$step;
  $ystep=$step;
  $zstep=$step;
}

unless($exact)
{
  # pad the xfm file
  $xstart-=$xstep;
  $ystart-=$ystep;
  $zstart-=$zstep;
  
  $xlen+=2*$xstep;
  $ylen+=2*$ystep;
  $zlen+=2*$zstep;
}

unless($exact)
{
  $xspace=POSIX::ceil($xlen/$xstep);
  $yspace=POSIX::ceil($ylen/$ystep);
  $zspace=POSIX::ceil($zlen/$zstep);
}

my $dimorder;

if($exact)
{
  $dimorder=join(',',split(' ',`mincinfo -vardims image $model`));
}

#disable _grid file compression!
delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};
$tmpdir= &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

if($in =~ /.mnc$|.mnc.gz$/)
{
  do_cmd('cp',$in,"$tmpdir/in_grid_0.mnc");
  open IN,">$tmpdir/in.xfm" or die;
  print IN <<XFM
MNI Transform File

Transform_Type = Grid_Transform;
Displacement_Volume = in_grid_0.mnc;
XFM
;
 close IN;
 $in="$tmpdir/in.xfm";
}

if($invert)
{
  do_cmd('xfminvert',$in,"$tmpdir/inverted.xfm");
  $in="$tmpdir/inverted.xfm";
}



my @arg=('xfm2def', '-xstep', $xstep,'-ystep', $ystep,'-zstep', $zstep,
  '-xstart', $xstart,'-ystart', $ystart,'-zstart', $zstart, 
  '-xnelements', $xspace,
  '-ynelements', $yspace, 
  '-znelements', $zspace,
  '-float', $in, $dimorder?"$tmpdir/grid.mnc":$out_grid,'-clobber');
do_cmd(@arg);

do_cmd('mincreshape','-dimorder',$dimorder,"$tmpdir/grid.mnc",$out_grid,'-clobber') if $dimorder;

open  OF,">$out" or die "Can't open ${out} for writing!\n";
print OF "MNI Transform File\n";
print OF "Transform_Type = Linear;\n";
print OF "Linear_Transform =\n";
print OF " 1 0 0 0\n 0 1 0 0\n 0 0 1 0;\n";
print OF "Transform_Type = Grid_Transform;\n";
print OF "Displacement_Volume = ";
print OF basename($out_grid);
print OF ";\n";
close OF;


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}



