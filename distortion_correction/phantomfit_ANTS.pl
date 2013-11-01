#! /usr/bin/env perl


############################# MNI Header #####################################
#@NAME       :  phantomfit.pl
#@DESCRIPTION:  estimate distortion field based on nonlinear registration using ANTs
#@COPYRIGHT  :
#              Vladimir S. Fonov  February, 2012
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
use warnings "all";
use Getopt::Tabular;
use File::Basename;
use File::Temp qw/ tempdir /;

# default minctracc parameters
my @def_ants_args = (
  '-t','SyN[0.25]',
  '--number-of-affine-iterations','0x0x0',
  '-i','100x100x20',
  # TODO: add regularization ?
   );

my($Help, $Usage, $me);
my(@opt_table, %opt, $outxfm, $outfile, @args, $tmpdir);

$me = &basename($0);
%opt = (
   'verbose'   => 1,
   'debug'     => 0,
   'clobber'   => 0,
   'fake'      => 0,
   'init_xfm'  => undef,
   'order'     => 5,
   'par'       => undef,
   'measure'   => undef,
   'step_iterations'=> undef,
   'min_step'  => 2,
   'limit'     => 0,
   'weight'    => 1,
   'keep'      => 1.0,
   'cyl'       => 0,
   'init'      => undef,
   'work_dir'  => undef,
   'pca'       => undef,
   'pcs'       => undef,
   );

$Help = <<HELP;
| $me does non-linear fitting between N files
| restricted by spherical harmonics functions
| 
| Problems or comments should be sent to: vladimir.fonov\@gmail.com
HELP

$Usage = "Usage: $me [options] source.mnc target.mnc fit_mask.mnc estimate_mask.mnc [source.mnc target.mnc source_mask.mnc target_mask.mnc...] output.xfm\n$me -help to list options\n\n";

@opt_table = (
   ["-verbose", "boolean", 0, \$opt{verbose},
      "be verbose" ],
   ["-debug", "boolean", 0, \$opt{debug},
      "for debugging" ],
   ["-clobber", "boolean", 0, \$opt{clobber},
      "clobber existing check files" ],
   ["-fake", "boolean", 0, \$opt{fake},
      "do a dry run, (echo cmds only)" ],
   ["-init_xfm", "string", 1, \$opt{init_xfm},
      "initial transformation (default identity)" ],
   ["-order","integer",1,\$opt{order},
      "Spherical harmonics order"],
   ["-par", "string", 1, \$opt{par},
      "Output parameters into a file" ],
   ["-measure", "string", 1, \$opt{measure},
      "Output measurements into a file" ], 
   ["-step_iterations", "integer", 1, \$opt{step_iterations},
      "Maximum number of iterations per step" ], 
   ["-min_step", "float", 1, \$opt{min_step},
      "Minimal step size for nonlinear registration (min 1.0)" ], 
   ["-limit", "boolean", 0, \$opt{limit},
      "Limit linear component to identity" ],
   ["-keep", "float", 1, \$opt{keep},
      "Fraction of points to keep for LTSQ algorithm (0-1]" ],
   ["-cylindric", "boolean", 0, \$opt{cyl},
      "Use cylindric functions instead of spherical" ],
   ["-init", "string", 1, \$opt{init},
      "Initial estimation" ],
   ["-work_dir", "string", 1, \$opt{work_dir},
      "Work directory (instead of temp), usefull for debugging" ],
   ["-pca", "string", 1, \$opt{pca},
      "Use pca rotation matrix" ],
   ["-pcs", "integer", 1, \$opt{pcs},
      "limit number of PCs" ],
   );

# Check arguments
&Getopt::Tabular::SetHelp($Help, $Usage);
&GetOptions (\@opt_table, \@ARGV) || exit 1;
die $Usage if  $#ARGV < 4;
my @source;
my @target;
my @fit_mask;
my @estimate_mask;
print "ARGV:",join(',',@ARGV),"\n";

my $minc_compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $minc_compress;

for(my $i=0;$i <= (($#ARGV-1)/4);$i++)
{
  push @source,$ARGV[$i*4];
  push @target,$ARGV[$i*4+1];
  push @fit_mask,$ARGV[$i*4+2];
  push @estimate_mask,$ARGV[$i*4+3];
}
print "Sources:",join(' ',@source),"\n";
print "Target:",join(' ',@target),"\n";
$outxfm=$ARGV[$#ARGV];
$outfile = $opt{output};

check_file($outxfm) unless $opt{clobber};
check_file($outfile) unless $opt{clobber} || !defined($outfile);
check_file($opt{par}) unless $opt{clobber} || !defined($opt{par});
#check_file($opt{measure}) unless $opt{clobber} || !defined($opt{measure});

# make tmpdir
unless($opt{work_dir})
{
  $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
} else {
  $tmpdir = $opt{work_dir};
  do_cmd('mkdir','-p',$tmpdir);
}

# set up filename base
my($i, $s_base, $t_base, $tmp_xfm, $tmp_grid, $tmp_source, $tmp_target, $prev_grid, $prev_xfm);
my @grids;

# a fitting we shall go...
for(my $k=0;$k<=$#source;$k++)
{
  $s_base = &basename($source[$k]);
  $s_base =~ s/\.gz$//;
  $s_base =~ s/\.mnc$//;

  $tmp_xfm  = "$tmpdir/$s_base\_${k}.xfm";
  $tmp_grid = "$tmpdir/$s_base\_${k}_grid_0.mnc";

  @args = ('mincANTS','3',@def_ants_args,
           '-m',"CC[$source[$k],$target[$k],1,2]",
           '-o',$tmp_xfm,
           '-x',$fit_mask[$k]);
            
  # transformation
  my $transform;
  
  &do_cmd(@args);
  push @grids,$tmp_grid;
}

regularize_grids(\@grids,\@estimate_mask,$opt{order},"$tmpdir/regularize.xfm",1);
cleanup_grids(\@grids) unless $opt{debug};

$prev_xfm = "$tmpdir/regularize.xfm";
$prev_grid = "$tmpdir/regularize_grid_0.mnc";
    
    
do_cmd('cp',"$tmpdir/regularize.par",$opt{par}) if $opt{par};
do_cmd('xfminvert',$prev_xfm,$outxfm);

sub do_cmd { 
   print STDOUT "@_\n" if $opt{verbose};
   if(!$opt{fake}){
      system(@_) == 0 or die "DIED: @_\n";
      }
   }
       
sub cleanup_xfms {
  my @xfms=@{$_[0]};
  for(my $i=0;$i<=$#xfms;$i++)
  {
    my $grid=$xfms[$i];
    $grid=~s/.xfm$/_grid_0.mnc/;
    do_cmd('rm','-f',$grid,$xfms[$i]);
  }
}

sub cleanup_grids {
  my @grids=@{$_[0]};
  for(my $i=0;$i<=$#grids;$i++)
  {
    do_cmd('rm','-f',$grids[$i]);
  }
}


sub regularize_grids {
  my ($_grid,$_mask,$order,$out,$invert,$out_par)=@_;
  my @grids=@{$_grid};
  my @masks=@{$_mask};
  die "XFMs and MASKs don't match!" if $#grids!=$#masks;
  
  my @args;
  if($opt{pca})
  {
    @args=('fit_harmonics_grids_regularize','--order',$opt{order},'--skip',2); # ANTS outputs grid @ 1mm
    push @args,'--cylindrical' if $opt{cyl};
    push @args,'--pca',$opt{pca};
    push @args,'--pcs',$opt{pcs} if $opt{pcs};
  } else {
    @args=($opt{cyl}?'c_fit_harmonics_grids':'fit_harmonics_grids','--order',$opt{order},'--skip',2);
    push(@args,"--limit") if $opt{limit};
    push(@args,"--keep",$opt{keep}) if $opt{keep};
    push(@args,'--iter',1) if $opt{keep}>0.99;
  }
  
  for(my $i=0;$i<=$#grids;$i++)
  {
    die "Can't find  grid file: $grids[$i]" unless -e $grids[$i];
    push (@args,$grids[$i],$masks[$i]);
  }
  
  push(@args,"$tmpdir/regularize.par",'--clobber');
  
  do_cmd(@args);
  @args=('par2xfm.pl',"$tmpdir/regularize.par",$out,'--clobber','--max',20,'--extent',400,'--step',4);
  push @args,'--noinvert'  if $invert;
  push @args,'--cylindric' if $opt{cyl};
  do_cmd(@args);
  do_cmd('cp',"$tmpdir/regularize.par",$out_par) if $out_par;
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}
