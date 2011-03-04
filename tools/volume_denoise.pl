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
my ($xy,$yz,$xz);
my $mask;
my $threads=1;
my $beta=1;

GetOptions(
	   'verbose' => \$verbose,
	   'clobber' => \$clobber,
     'threads=n' => \$threads,
     'beta=f' => \$beta
     );
     
my $help=<<HELP
Volume denoising script based on 
Pierrick Coupe, Jose V. Manjon, Elias Gedamu, Douglas L. Arnold,
Montserrat Robles, D. Louis Collins: An Object-Based Method for Rician
Noise Estimation in MR Images. MICCAI (1) 2009: 601-608

and

P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.
An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic
Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441,
Avril 2008
.
Usage: $me <input> <output> [--verbose --clobber --threads <n, default 1> --beta <n, 0<beta<=1, default 1>]\n
HELP
;

die $help if $#ARGV<1;

my ($in,$out)=@ARGV;

check_file($out) if !$clobber;

#1 estimate noise level 

my $noise=`noise_estimate $in`;
chomp($noise);

print "Noise level=$noise\n" if $verbose;

# check if the volume have more or less uniform step size

my @steps=split(/\n/,`mincinfo -attvalue xspace:step -attvalue yspace:step  -attvalue zspace:step $in`);

my $i;
for($i=0;$i<=$#steps;$i+=1)
{
  $steps[$i]=abs($steps[$i]);
}
my @steps_ = sort {$a <=> $b} @steps;
my $aniso=abs($steps_[0]/$steps_[2])<0.5?1:0; # nonuniform step size?

if($aniso) # add -beta $beta
{
  print "Volume has anisotropic resolution\n";
  do_cmd('mincnlm',$in,$out,'-clobber','-sigma',$noise,'-aniso','-mt',$threads,'-v ','2','-d','9','-w','2');
} else {
  do_cmd('mincnlm',$in,$out,'-clobber','-sigma',$noise,'-mt',$threads,'-w','2');
}

sub do_cmd { 
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
      system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}