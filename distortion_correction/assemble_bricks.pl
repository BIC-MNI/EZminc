#!/usr/bin/env perl

############################# MNI Header #####################################
#@NAME       :  assemble_bricks.pl
#@DESCRIPTION:  create lego phantom image based on location of bricks
#@COPYRIGHT  :
#              Vladimir S. Fonov  April 2012, April 2018
#              Montreal Neurological Institute, McGill University.
#              Permission to use, copy, modify, and distribute this
#              software and its documentation for any purpose and without
#              fee is hereby granted, provided that the above copyright
#              notice appear in all copies.  The author and McGill University
#              make no representations about the suitability of this
#              software for any purpose.  It is provided "as is" without
#              express or implied warranty.
###############################################################################

use strict;
use File::Basename;
use File::Temp qw/tempdir/;
use Getopt::Long;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $brick='';
my $me = &basename($0);
my $out;
my $like;
my $mask;
my $cont=0;
my $tag;
my $tmpdir;

GetOptions (    
	        "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "brick=s"   => \$brick,
          "like=s"    => \$like,
          "mask=s"    => \$mask,
          "tmp=s"     => \$tmpdir
          );
          
die "Usage: $me <bricks.tag> <output_model_ideal.mnc> <output_model_bricks.mnc> --verbose --clobber --brick <brick> --mask <brick_mask> --like <sample> \n"  if $#ARGV<0;

my ($tag, $out_phantom, $out_bricks)=@ARGV;

die "Need a brick and mask!\n" unless $brick || $mask;
die "Specify Output please\n" unless $out_phantom && $out_bricks;
die "Specify an example (--like)\n" unless $like;

check_file($out) unless $clobber;
$tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 ) unless $tmpdir;

my @xfm;

my @tags=parse_tags($tag);
do_cmd('param2xfm', '-rotations',0 ,0, 90, "$tmpdir/rotate90.xfm");
my $i;
for $i(0..$#tags)
{
	my ($x,$y,$z,$t)=@{$tags[$i]};
	my $out_tmp="$tmpdir/$i.xfm";
	my $out_xfm="$tmpdir/_${t}_$i.xfm";
	do_cmd('param2xfm','-translation', $x, $y, $z,$out_tmp,'-clobber');
	if($t==1)
	{
		do_cmd('xfmconcat',"$tmpdir/rotate90.xfm",$out_tmp,$out_xfm);
	} else {
		do_cmd('cp', $out_tmp, $out_xfm);
	}
	push (@xfm, $out_xfm);
}


my $comp=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $comp;

my $b;
my $j=1;
my @bricks;
do_cmd('minccalc','-expression', 0, $like, "$tmpdir/brick_tmp.mnc", '-byte', '-labels');

# 1 create continious phantom  
foreach $b(@xfm)
{
  do_cmd('mincresample', $brick, "$tmpdir/brick_${j}.mnc", '-transform', $b, '-like', $like,'-clobber');
  do_cmd('mincmath','-max', "$tmpdir/brick_${j}.mnc", "$tmpdir/brick_tmp.mnc", "$tmpdir/brick_tmp2.mnc", '-clobber');
  do_cmd('mv', "$tmpdir/brick_tmp2.mnc", "$tmpdir/brick_tmp.mnc");
  do_cmd('rm', '-f', "$tmpdir/brick_${j}.mnc");
  $j++;
}
do_cmd('mincreshape', "$tmpdir/brick_tmp.mnc", $out_phantom);
do_cmd('rm', '-f', "$tmpdir/brick_tmp.mnc");
do_cmd('minccalc', '-expression', 0, $like, "$tmpdir/brick_tmp.mnc", '-byte', '-labels');


unless($mask) {
    $mask="$tmpdir/brick.mnc";
    do_cmd('itk_morph', '--bimodal', $brick, $mask);
}

# 2 create bricks view
$j=1;
foreach $b(@xfm)
{
  do_cmd('mincresample', $mask, "$tmpdir/brick_t.mnc", '-nearest', '-transform', $b, '-like', $like, '-clobber', '-labels');
  do_cmd('minccalc', '-expression',"A[0]*$j", "$tmpdir/brick_t.mnc","$tmpdir/brick_${j}.mnc",'-labels');
  do_cmd('mincmath', '-max', "$tmpdir/brick_${j}.mnc", "$tmpdir/brick_tmp.mnc", "$tmpdir/brick_tmp2.mnc", '-clobber','-labels');
  do_cmd('mv', "$tmpdir/brick_tmp2.mnc","$tmpdir/brick_tmp.mnc");
  do_cmd('rm', '-f', "$tmpdir/brick_${j}.mnc");
  $j++;
}

$ENV{MINC_COMPRESS}=$comp if $comp;

do_cmd('mincreshape',"$tmpdir/brick_tmp.mnc",$out_bricks, '-image_range',0,$j, '-valid_range',0,$j);


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}
sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}


#read tags from the file, with additional type parameter
sub parse_tags {
  my ($tag)=@_;
  open TAG,"<$tag" or die "Can't open $tag\n";
  my $line;
  my @tags;
  my $started=0;
  foreach  $line(<TAG>)
  {
    if(!$started)
    {
      if($line =~ /.*Points =/) 
      {
          $started=1;
          next;
      }
      next;
    }
    if($line =~ /.*;/) # this is the last line
    {
      $line =~ s/;//g;
    }
    my @c=split(/\s/, $line);
    shift @c unless $c[0]; #protection against empty first parameter
    if($#c==3) {
      my $l=$c[3];
      $l=~ s/\"//g;
      push(@tags, [$c[0], $c[1], $c[2], $l+0]);
    } else {
      my $l=$c[3];
      $l=~ s/\"//g;
      push(@tags, [$c[0], $c[1], $c[2], $l+0]);
    }
  }

  close TAG;
  return @tags;
}

