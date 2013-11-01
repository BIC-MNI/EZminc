#! /usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use POSIX qw(floor);

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $resample='tricubic';
my $step=1.0;
my $transform;
my $datatype;
my $invert_transform=0;
my $debug=0;

GetOptions (    
          "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "step=f"    => \$step,
          "resample=s" =>\$resample,
          "transform=s" =>\$transform,
          "datatype=s" =>\$datatype,
          'invert_transform' => \$invert_transform
          );
          
die "Usage: $me <scan_in> <scan_out> \n[--verbose\n --clobber\n --step <f>\n --resample <type>\n --transform <transform.xfm>\n --datatype <type>\n --invert_transform ]\n" if $#ARGV<1;

my ($in,$out)=@ARGV;

check_file($out) unless $clobber;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

my ($xspace,$yspace,$zspace, $xstart,$ystart,$zstart, $xstep,$ystep,$zstep)= split(/\n/, 
 `mincinfo  -dimlength xspace -dimlength yspace -dimlength zspace -attvalue xspace:start -attvalue yspace:start -attvalue zspace:start -attvalue xspace:step -attvalue yspace:step -attvalue zspace:step $in`);

my @dircos=split(/\n/,`mincinfo -attvalue xspace:direction_cosines -attvalue yspace:direction_cosines -attvalue zspace:direction_cosines $in`);

my @xdircos;
my @ydircos;
my @zdircos;

if($#dircos<2)
{
  @xdircos=(1,0,0);
  @ydircos=(0,1,0);
  @zdircos=(0,0,1);
} else {
  @xdircos=split(/\s/,$dircos[0]);
  @ydircos=split(/\s/,$dircos[1]);
  @zdircos=split(/\s/,$dircos[2]);
}

#I hate perl!
$xstart=$xstart*1.0;
$ystart=$ystart*1.0;
$zstart=$zstart*1.0;

$xstep=$xstep*1.0;
$ystep=$ystep*1.0;
$zstep=$zstep*1.0;


if($xstep<0)
{
  $xstart+=$xstep*($xspace-1);
  $xstep=-$xstep;
}

if($ystep<0)
{
 $ystart+=$ystep*($yspace-1);
 $ystep=-$ystep;
}

if($zstep<0)
{
  $zstart+=$zstep*($zspace-1); 
  $zstep=-$zstep;
}

if($debug)
{
  print "$xstart,$ystart,$zstart\n";
  print "$xspace,$yspace,$zspace\n";
  print join(',',@xdircos),"\n";
  print join(',',@ydircos),"\n";
  print join(',',@zdircos),"\n";
}

my @mat;
my $i;
for ($i=0;$i<3;$i+=1)
{
  push(@mat,[ $xdircos[$i],$ydircos[$i],$zdircos[$i]]);
}

my @edge;
push(@edge,mult_add(\@mat, [0,0,0],                  [$xstart,$ystart,$zstart]));
push(@edge,mult_add(\@mat, [$xspace*$xstep,0,0],            [$xstart,$ystart,$zstart]));
push(@edge,mult_add(\@mat, [0,$yspace*$ystep,0],            [$xstart,$ystart,$zstart]));
push(@edge,mult_add(\@mat, [0,0,$zspace*$zstep],            [$xstart,$ystart,$zstart]));
push(@edge,mult_add(\@mat, [$xspace*$xstep,0,$zspace*$zstep],      [$xstart,$ystart,$zstart]));
push(@edge,mult_add(\@mat, [$xspace*$xstep,$yspace*$ystep,0],      [$xstart,$ystart,$zstart]));
push(@edge,mult_add(\@mat, [0,$yspace*$ystep,$zspace*$zstep],      [$xstart,$ystart,$zstart]));
push(@edge,mult_add(\@mat, [$xspace*$xstep,$yspace*$ystep,$zspace*$zstep],[$xstart,$ystart,$zstart]));

my (@xedge,@yedge,@zedge);

foreach ($i=0;$i<=$#edge;$i+=1)
{
  push(@xedge,$edge[$i][0]);
  push(@yedge,$edge[$i][1]);
  push(@zedge,$edge[$i][2]);
}

if($debug)
{
  print join(',',@xedge),"\n";
  print join(',',@yedge),"\n";
  print join(',',@zedge),"\n\n";
}

my @xxedge=sort {$a<=>$b} @xedge;
my @yyedge=sort {$a<=>$b} @yedge;
my @zzedge=sort {$a<=>$b} @zedge;

#print join(',',@xxedge),"\n";
#print join(',',@yyedge),"\n";
#print join(',',@zzedge),"\n";
$xstart=$xxedge[0];
$ystart=$yyedge[0];
$zstart=$zzedge[0];

my $xlen=$xxedge[$#xxedge]-$xstart;
my $ylen=$yyedge[$#yyedge]-$ystart;
my $zlen=$zzedge[$#zzedge]-$zstart;

my @arg=('mincresample', '-'.$resample, $in, $tmpdir."/resampled.mnc", 
    '-dircos', 1, 0, 0, 0, 1, 0, 0, 0 ,1 , 
    '-step', $step, $step, $step, 
    '-start', $xstart, $ystart, $zstart, 
    '-nelements', POSIX::floor($xlen/$step), POSIX::floor($ylen/$step), POSIX::floor($zlen/$step));

push @arg,('-transform',$transform,'-tfm_input_sampling') if $transform ; 
push @arg,('-invert_transformation') if $invert_transform;
push @arg,'-'.$datatype if $datatype;

my $compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $compress;

do_cmd(@arg);

$ENV{MINC_COMPRESS}=$compress if $compress;;
do_cmd('mincreshape','-dimorder','zspace,yspace,xspace',$tmpdir."/resampled.mnc", $out,'-clobber');

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}

sub mult {
  my ($mat_,$vec_)=@_;
  my @mat=@{$mat_};
  my @vec=@{$vec_};
  my ($i,$j);
  my @res=(0,0,0);
  for($i=0;$i<3;$i+=1){
    for($j=0;$j<3;$j+=1) {
      $res[$i]+=$mat[$i][$j] * $vec[$j];
    }
  }
  return \@res;
}

sub mult_add {
  my ($mat_,$vec_,$add_)=@_;
  my @mat=@{$mat_};
  my @vec=@{$vec_};
  my @add=@{$add_};
  my ($i,$j);
  my @res=(0.0,0.0,0.0);
  for($i=0;$i<3;$i+=1){
    for($j=0;$j<3;$j+=1) {
      $res[$i]+=$mat[$i][$j] * ($vec[$j]+$add[$j]) ;
    }
    #$res[$i]+= $add[$i];
  }
  
  return \@res;
}
