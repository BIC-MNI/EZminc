#! /usr/bin/env perl
#

use strict;
use warnings "all";
use Getopt::Tabular;
use File::Basename;
use File::Temp qw/ tempdir /;

# default minctracc parameters
my @def_minctracc_args = (
#   '-debug',
   '-clobber',
   '-nonlinear', 'corrcoeff',
   '-sub_lattice', 6
   );

my @conf = (

   {'step'         => 4,
    'blur_fwhm'    => 8,
    'iterations'   => 6,
    },

    {'step'         => 4,
     'blur_fwhm'    => 4,
     'iterations'   => 5,
    },

    {'step'         => 3,
     'blur_fwhm'    => 4,
     'iterations'   => 3,
    },
    
    {'step'         => 2,
     'blur_fwhm'    => 2,
     'iterations'   => 3,
    },
    
    {'step'         => 1,
     'blur_fwhm'    => 2,
     'iterations'   => 3,
    },
    
   );

my($Help, $Usage, $me);
my(@opt_table, %opt, $outxfm, $outfile, @args, $tmpdir);

$me = &basename($0);
%opt = (
   'verbose'   => 0,
   'debug'     => 0,
   'clobber'   => 0,
   'fake'      => 0,
   'init_xfm'  => undef,
   'order'     => 3,
   'par'       => undef,
   'measure'   => undef,
   'step_iterations'=> undef,
   'min_step'  => 2,
   'limit'     => 0,
   'weight'    => 1,
   'stiffness' => 1,
   'similarity'=> 0.3,
   'keep'      => 1.0,
   'cyl'       => 0,
   'init'      => undef,
   'work_dir'  => undef
   );

$Help = <<HELP;
| $me does hierachial non-linear fitting between N files
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
   ["-weight", "float", 1, \$opt{weight},
      "Wight parameter of minctracc" ],
   ["-stiffness", "float", 1, \$opt{stiffness},
      "stiffness parameter of minctracc" ],
   ["-similarity", "float", 1, \$opt{similarity},
      "similarity parameter of minctracc" ],
   ["-keep", "float", 1, \$opt{keep},
      "Fractioon of points to keep for LTSQ algorithm (0-1]" ],
   ["-cylindric", "boolean", 0, \$opt{cyl},
      "Use cylindric functions instead of spherical" ],
   ["-init", "string", 1, \$opt{init},
      "Initial estimation" ],
   ["-work_dir", "string", 1, \$opt{work_dir},
      "Work directory (instead of temp), usefull for debugging" ],
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
check_file($opt{measure}) unless $opt{clobber} || !defined($opt{measure});

if($opt{measure})
{
  open MEASURE,">$opt{measure}" ;
  print MEASURE "step,blur,measure\n";
}
# make tmpdir
unless($opt{work_dir})
{
  $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );
} else {
  $tmpdir = $opt{work_dir};
  do_cmd('mkdir','-p',$tmpdir);
}

if($opt{init})
{
  do_cmd('xfminvert',$opt{init},"$tmpdir/init.xfm");
}


# set up filename base
my($i, $s_base, $t_base, $tmp_xfm, $tmp_source, $tmp_target, $prev_xfm);

#if(defined $opt{init_xfm})
#{
# regularize_xfms([$opt{init_xfm}],[$opt{source_mask}],$opt{order},"$tmpdir/regularize.xfm",1);
#}

# a fitting we shall go...
for ($i=0; $i<=$#conf; $i++) {
  last if $conf[$i]{step}<$opt{min_step};
  
  for(my $k=0;$k<=$#source;$k++)
  {
    $s_base = &basename($source[$k]);
    $s_base =~ s/\.gz$//;
    $s_base =~ s/\.mnc$//;
    $t_base = &basename($target[$k]);
    $t_base =~ s/\.gz$//;
    $t_base =~ s/\.mnc$//;

    # set up intermediate files
    $tmp_source = "$tmpdir/S_$s_base\_$conf[$i]{blur_fwhm}";
    $tmp_target = "$tmpdir/T_$t_base\_$conf[$i]{blur_fwhm}";
    $conf[$i]{iterations}=$opt{step_iterations} if $opt{step_iterations};
    print STDOUT "-+-[$i]\n".
                  " | step:           $conf[$i]{step}\n".
                  " | blur_fwhm:      $conf[$i]{blur_fwhm}\n".
                  " | iterations:     $conf[$i]{iterations}\n".
                  " | source:         $tmp_source\n".
                  " | target:         $tmp_target\n".
                  "\n";
   
   # blur the source and target files if required
     if(!-e "$tmp_source\_blur.mnc"){
        &do_cmd('mincblur','-no_apodize', '-fwhm', $conf[$i]{blur_fwhm}, $source[$k], $tmp_source);
        }
     if(!-e "$tmp_target\_blur.mnc"){
        &do_cmd('mincblur','-no_apodize', '-fwhm', $conf[$i]{blur_fwhm}, $target[$k], $tmp_target);
        }
     my $j;
   }
   
  for(my $j=0;$j<$conf[$i]{iterations};$j++)
  {
     my @xfms;
     
    if($opt{measure}) 
    {
      my @mes;
      if($i == 0 && $j == 0 )
      {
        @mes=split(/Distance\:/,`nl_distance.pl $source[0] $target[0] --mask $fit_mask[0]`);
      } else {
        # now we need to resample the source & mask
        do_cmd('mincresample',$source[0],'-like',$target[0],"$tmpdir/measure.mnc",'-transform',$prev_xfm);
        do_cmd('mincresample',$fit_mask[0],'-like',$target[0],"$tmpdir/measure_mask.mnc",'-transform',$prev_xfm,'-nearest');
        @mes=split(/Distance\:/,`nl_distance.pl $tmpdir/measure.mnc $target[0] --mask $tmpdir/measure_mask.mnc`);
        do_cmd('rm','-f',"$tmpdir/measure.mnc","$tmpdir/measure_mask.mnc");
      }
      my $m=pop @mes;
      chomp($m);
      print MEASURE "$conf[$i]{step},$conf[$i]{blur_fwhm},$m\n"
    }
     
    for(my $k=0;$k<=$#source;$k++)
    {
      $s_base = &basename($source[$k]);
      $s_base =~ s/\.gz$//;
      $s_base =~ s/\.mnc$//;
      $t_base = &basename($target[$k]);
      $t_base =~ s/\.gz$//;
      $t_base =~ s/\.mnc$//;
      $tmp_source = "$tmpdir/S_$s_base\_$conf[$i]{blur_fwhm}";
      $tmp_target = "$tmpdir/T_$t_base\_$conf[$i]{blur_fwhm}";

      $tmp_xfm = "$tmpdir/$s_base\_${i}_${j}_${k}.xfm";
      # set up registration
      @args = ('minctracc',  @def_minctracc_args,
            '-iterations', 2,
            '-step', $conf[$i]{step}, $conf[$i]{step}, $conf[$i]{step},
            '-lattice_diam', $conf[$i]{step} * 3, 
                             $conf[$i]{step} * 3, 
                             $conf[$i]{step} * 3,
            '-similarity_cost_ratio',$opt{similarity},
            '-weight',$opt{weight},
            '-stiffness',$opt{stiffness});

      # transformation
      my $transform;
      if($i == 0 && $j == 0 ){ # no init xfm for now!
        if($opt{init})
        {
          push(@args, '-transformation',"$tmpdir/init.xfm" );
        } else {
          push(@args, '-identity');
        }
      } else {
        push(@args, '-transformation', $prev_xfm );
      }

      push(@args,'-no_super') if $conf[$i]{step}<4;

      # masks
      push(@args, '-source_mask',$fit_mask[$k]) ;
      push(@args, '-model_mask', $fit_mask[$k]);

      # add files and run registration
      push(@args, "$tmp_source\_blur.mnc", "$tmp_target\_blur.mnc", $tmp_xfm);
      &do_cmd(@args);
      push @xfms,$tmp_xfm;
    }
    regularize_xfms(\@xfms,\@estimate_mask,$opt{order},"$tmpdir/regularize_${i}_${j}.xfm",1);
    cleanup_xfms(\@xfms) unless $opt{debug};
    $prev_xfm = "$tmpdir/regularize_${i}_${j}.xfm";
  }
}
#regularize_xfm($prev_xfm,$opt{source_mask},$opt{order},$outxfm,0,$opt{par});
do_cmd('cp',"$tmpdir/regularize.par",$opt{par}) if $opt{par};
do_cmd('xfminvert',$prev_xfm,$outxfm);
close MEASURE if $opt{measure};

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

sub regularize_xfms{
  my ($_xfm,$_mask,$order,$out,$invert,$out_par)=@_;
  my @xfms=@{$_xfm};
  my @masks=@{$_mask};
  die "XFMs and MASKs don't match!" if $#xfms!=$#masks;
  
  my @args=($opt{cyl}?'c_fit_harmonics_grids':'fit_harmonics_grids','--order',$opt{order});
  for(my $i=0;$i<=$#xfms;$i++)
  {
    my $grid=$xfms[$i];
    $grid=~s/.xfm$/_grid_0.mnc/;
    die "Can't find  grid file: $grid" unless -e $grid;
    push (@args,$grid,$masks[$i]);
  }
  push(@args,"$tmpdir/regularize.par",'--iter',10,'--clobber');
  push(@args,"--limit") if $opt{limit};
  push(@args,"--keep",$opt{keep}) if $opt{keep};
  
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