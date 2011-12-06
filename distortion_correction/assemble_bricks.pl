#!/usr/bin/env perl
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
          "out=s"     => \$out,
          "like=s"    => \$like,
          "mask=s"    => \$mask,
          "cont"      => \$cont,
          "tmp=s"     => \$tmpdir
          );
          
die "Usage: $me <bricks.tag> <output_model.mnc> --verbose --clobber --brick <brick> --out <out> --mask <brick_mask> --like <sample> --cont \n"  if $#ARGV<0;

my ($tag,$out)=@ARGV;

die "Need a brick or mask!\n" unless $brick || $mask;
die "Specify Output please\n" unless $out;
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
		do_cmd('cp', $out_tmp,$out_xfm);
	}
	push (@xfm,$out_xfm);
}

my $b;
my $j=1;
my @bricks;
do_cmd('minccalc','-expression',0,$like,"$tmpdir/brick_tmp.mnc",'-byte');
if($cont) {
  
  foreach $b(@xfm)
  {
    do_cmd('mincresample',$brick,"$tmpdir/brick_${j}.mnc",'-transform',$b,'-like', $like,'-clobber');
    do_cmd('mincmath','-max',"$tmpdir/brick_${j}.mnc","$tmpdir/brick_tmp.mnc","$tmpdir/brick_tmp2.mnc",'-clobber');
    do_cmd('mv',"$tmpdir/brick_tmp2.mnc","$tmpdir/brick_tmp.mnc");
    do_cmd('rm','-f',"$tmpdir/brick_${j}.mnc");
    $j++;
  }
} else {
  unless($mask) {
      $mask="$tmpdir/brick.mnc";
      do_cmd('itk_morph', '--bimodal', $brick, $mask);
  }
  foreach $b(@xfm)
  {
    do_cmd('mincresample',$mask,"$tmpdir/brick_t.mnc",'-nearest','-transform',$b,'-like', $like,'-clobber');
    do_cmd('minccalc','-expression',"A[0]*$j","$tmpdir/brick_t.mnc","$tmpdir/brick_${j}.mnc");
    do_cmd('mincmath','-max',"$tmpdir/brick_${j}.mnc","$tmpdir/brick_tmp.mnc","$tmpdir/brick_tmp2.mnc",'-clobber');
    do_cmd('mv',"$tmpdir/brick_tmp2.mnc","$tmpdir/brick_tmp.mnc");
    do_cmd('rm','-f',"$tmpdir/brick_${j}.mnc");
    $j++;
  }
}

#do_cmd('mincmath','-max',@bricks,$out,'-clobber');
do_cmd('cp',"$tmpdir/brick_tmp.mnc",$out);


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

