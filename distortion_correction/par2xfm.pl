#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $me=basename($0);
my $spacing=2;
my $max=5.0;
my $noinvert=0;
my $extent=300;
my $history=localtime() .">>> ".$me." ".join(' ',@ARGV);
my $cylindric=0;

GetOptions (    
	  "verbose"   => \$verbose,
          "clobber"   => \$clobber,
          "spacing=f" => \$spacing,
          "max=f"     => \$max,
          "noinvert"  => \$noinvert,
          "extent=f"  => \$extent,
	  "step=f"    => \$spacing,
	  "cylindric" => \$cylindric); 
          
die "Programm usage: $me <par_in> <xfm_out> [--clobber] [--verbose]  [--spacing <n>] [--max <f>] [--noinvert] [--extent <mm>] [--cylindric]\n" if $#ARGV<1;

my ($par,$xfm)=@ARGV;

my $output_grid=$xfm;
$output_grid=~s/\.xfm//;
$output_grid.='_grid_0.mnc';

check_file($output_grid) if !$clobber;
check_file($xfm) if !$clobber;

do_cmd($cylindric?'c_param2grid':'param2grid',$par,$output_grid,'--spacing', $spacing,'--clobber', '--max', $max,'--extent',$extent);

#5 make xfm file
open  OF,">$xfm" or die "Can't open ${xfm} for writing!\n";
print OF "MNI Transform File\n";
print OF "Transform_Type = Linear;\nLinear_Transform =\n 1 0 0 0\n 0 1 0 0\n 0 0 1 0;\n";
print OF "Transform_Type = Grid_Transform;\n";
print OF "Invert_Flag = True;\n" unless $noinvert;
print OF "Displacement_Volume =";
print OF basename($output_grid);
print OF ";";
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
