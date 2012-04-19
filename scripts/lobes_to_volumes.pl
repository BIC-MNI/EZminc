#!/usr/bin/env perl

##############################################################################
#
# Vladimir S. Fonov July 2007
# McConnell Brain Imaging Centre, 
# Montreal Neurological Institute, 
# McGill University
##############################################################################

use strict;
use Getopt::Long;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;


my $me= basename($0);
my $verbose     = 1;
my $xfm;
my $infile_lobes;

GetOptions(
	   'xfm=s'   => \$xfm
	   );

if($#ARGV < 0){ die "Usage: $me <input lobes file> [--xfm <xfm file>]\n"; }

# Get the arguments.
$infile_lobes = $ARGV[0];

my @types = ("total", "wm" , "gm", "csf");
my @lobes = ("FrontalLobe_Right", "ParietalLobe_Right", "TemporalLobe_Right", "OccipitalLobe_Right", "FrontalLobe_Left", "ParietalLobe_Left", "TemporalLobe_Left", "OccipitalLobe_Left");

my %lobe_map = (30=>'frontal_left_wm',
				210=>'frontal_left_gm',
				17=>'frontal_right_wm',
				211=>'frontal_right_gm',
				83=>'temporal_left_wm',
				218=>'temporal_left_gm',
				59=>'temporal_right_wm',
				219=>'temporal_right_gm',
				57=>'parietal_left_wm',
				6=>'parietal_left_gm',
				105=>'parietal_right_wm',
				2=>'parietal_right_gm',
				73=>'occipital_left_wm',
				8=>'occipital_left_gm',
				45=>'occipital_right_wm',
				4=>'occipital_right_gm',
				67=>'cerebellum_left',
				76=>'cerebellum_right',
				20=>'brainstem',
				3=>'lateral_ventricle_left',
				9=>'lateral_ventricle_right',
				232=>'3rd_ventricle',
				233=>'4th_ventricle',
				255=>'extracerebral_CSF',
				39=>'caudate_left',
				53=>'caudate_right',
				14=>'putamen_left',
				16=>'putamen_right',
				102=>'thalamus_left',
				203=>'thalamus_right',
				33=>'subthalamic_nucleus_left',
				23=>'subthalamic_nucleus_right',
				12=>'globus_pallidus_left',
				11=>'globus_pallidus_right',
				29=>'fornix_left',
				254=>'fornix_right',
				28=>'skull');


# Get and print the scale factor.
#my $scalefactor = GetScaleFactor($infile_xfm);
my $scalefactor = 1.0;
$scalefactor = GetScaleFactor($xfm) if $xfm;
print "scale ${scalefactor}\n";
my @results = split(/\n/, `print_all_labels ${infile_lobes}`);
my $line;
foreach $line(@results) {
	chomp $line;
	my ($dummy,$label, $value) = split(/\s/, $line);
	$value = $value / $scalefactor;
	next unless $lobe_map{$label};
#	print "Lobe: ${lobe_map{$label}} = ${value}\n";
	print "${lobe_map{$label}} ${value}\n";
}

sub GetScaleFactor
{
  my $in_xfmfile = $_[0];
  
  my $scale_factor = 1.0;
        
  if (! -e $in_xfmfile ) {
     warn "Missing file $in_xfmfile\n";
  }
  else {
    my $scale = `xfm2param $in_xfmfile | grep scale`;
    my ($d, $scale_x, $scale_y, $scale_z) = split(/\s+/, $scale);
    $scale_factor = $scale_x * $scale_y * $scale_z;
 }
                                                      
 return $scale_factor;
}
                                                            
