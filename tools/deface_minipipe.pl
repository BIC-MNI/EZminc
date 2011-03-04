#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;

my $fake=0;
my $verbose=0;
my $clobber=0;
my $amp=1.0;
my $fwhm=10;
my $me=basename($0);
my $keep_tmp=0;
my $mask;
my $byte;
my $edge_smooth;
my $dual_echo=0;
my $model_dir;
my $model_name;
my $nonlinear=0;

my $amp=12;
my $fwhm=6;

my $no_int_norm=0;
my $watermark=0;
my $t1w_xfm;
my $t2w_xfm;
my $pdw_xfm;
my $brain_mask;
my $keep_real_range;

GetOptions (    
  "verbose"       => \$verbose,
  "clobber"       => \$clobber,
  "dual-echo"     => \$dual_echo,
  "model-dir=s"   => \$model_dir,
  "model=s"       => \$model_name,
  "nonlinear"     => \$nonlinear,
  'no-int-norm'   => \$no_int_norm,
  'watermark'     => \$watermark,
  't1w-xfm=s'     => \$t1w_xfm,
  't2w-xfm=s'     => \$t2w_xfm,
  'pdw-xfm=s'     => \$pdw_xfm,
  'brain-mask=s'  => \$brain_mask,
  'keep-tmp'      => \$keep_tmp,
  'keep-real-range' => \$keep_real_range
  
); 

die <<HELP
Usage: $me <T1w> [T2w] [PDw] <output_base>
 --model-dir <model directory>
 --model     <model name>
 [--verbose
  --clobber
  --dual-echo - assume dual echo T2/PD
  --nonlinear  perform quick nonlinear registration to compensate for different head shape (i.e for small kids)
  --watermark watermark output files
  --no-int-norm - don't normalize output file to 0-4095 range, usefull only for watermarking
  --t1w-xfm <t1w.xfm>
  --t2w-xfm <t2w.xfm>
  --pdw-xfm <pdw.xfm>
  --keep-real-range - keep the real range of the data the same 
  ]
HELP
if $#ARGV<1;

my $output_base=pop @ARGV;

my ($t1w,$t2w,$pdw)=@ARGV;

my $out_grid="${output_base}_deface_grid_0.mnc";

my $out_t1w="${output_base}_t1w.mnc";
my $out_t2w="${output_base}_t2w.mnc";
my $out_pdw="${output_base}_pdw.mnc";

#check_file($out_grid) unless $clobber;

check_file($out_t1w) unless $clobber;
check_file($out_t2w) if !$clobber||$t2w ;
check_file($out_pdw) if !$clobber||$pdw ;

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => !$keep_tmp );

my $compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};

my $model_t1w="$model_dir/$model_name.mnc";
my $model_t2w=$model_t1w;
my $model_pdw=$model_t1w;
my $model_face="$model_dir/${model_name}_face_mask.mnc";

$model_t2w=~s/t1/t2/;
$model_pdw=~s/t1/pd/;

my $model_mask="$model_dir/${model_name}_mask.mnc";
#fix irregular sampling
$t1w=fix_sampling($t1w);
$t2w=fix_sampling($t2w) if $t2w;
$pdw=fix_sampling($pdw) if $pdw;

# perform NU correct & clamp
correct($t1w,$model_t1w,"$tmpdir/clp_t1w.mnc") unless $t1w_xfm && $brain_mask;
correct($t2w,$model_t2w,"$tmpdir/clp_t2w.mnc") if $t2w && !($t1w_xfm&&$brain_mask);
correct($pdw,$model_pdw,"$tmpdir/clp_pdw.mnc") if $pdw && !($t1w_xfm&&$brain_mask);

#stx registration
unless($t1w_xfm&&$brain_mask)
{
  unless($t1w_xfm)
  {
    do_cmd('bestlinreg.pl','-lsq9',"$tmpdir/clp_t1w.mnc",$model_t1w,"$tmpdir/t1w_stx.xfm") ;
    $t1w_xfm="$tmpdir/t1w_stx.xfm";
  } 
  #T2 to T1 registration

  if($t2w && ! $t2w_xfm)
  {
      do_cmd('mritoself','-mi','-lsq6','-close','-nothreshold',
          "$tmpdir/clp_t2w.mnc","$tmpdir/clp_t1w.mnc",
          "$tmpdir/t2t1.xfm");
          
      do_cmd('xfmconcat',"$tmpdir/t2t1.xfm","$tmpdir/t1w_stx.xfm","$tmpdir/t2w_stx.xfm");
      $t2w_xfm="$tmpdir/t2w_stx.xfm";
  }

  if($pdw && !$pdw_xfm)
  {
    if($dual_echo)
    {
      do_cmd('cp',"$tmpdir/t2w_stx.xfm","$tmpdir/pdw_stx.xfm");
    } else {
      do_cmd('mritoself','-mi','-lsq6','-close','-nothreshold',
          "$tmpdir/clp_pdw.mnc","$tmpdir/clp_t1w.mnc",
          "$tmpdir/pdt1.xfm");
      do_cmd('xfmconcat',"$tmpdir/pdt1.xfm","$tmpdir/t1w_stx.xfm","$tmpdir/pdw_stx.xfm");    
    }
    $pdw_xfm="$tmpdir/pdw_stx.xfm";
  }
}

unless($brain_mask)
{
  do_cmd('itk_resample',"$tmpdir/clp_t1w.mnc","$tmpdir/stx_t1w.mnc",'--transform',$t1w_xfm,'--like',$model_t1w);
  do_cmd('mincbet',"$tmpdir/stx_t1w.mnc","$tmpdir/stx_brain",'-m','-n');
  $brain_mask="$tmpdir/stx_brain_mask.mnc";
}

if( -e $out_grid && $clobber)
{
  do_cmd('rm','-f',$out_grid);
}

unless( -e $out_grid )
{
  if($nonlinear)
  {
    do_cmd('itk_resample',"$tmpdir/clp_t1w.mnc","$tmpdir/stx_t1w.mnc",'--transform',$t1w_xfm,'--like',$model_t1w) if( ! -e "$tmpdir/stx_t1w.mnc");
    do_cmd('nlfit_s','-level',8,"$tmpdir/stx_t1w.mnc",$model_t1w,"$tmpdir/nl.xfm");

    do_cmd('mincresample','-nearest',$model_face,'-transform',"$tmpdir/nl.xfm",'-use_input_sampling',"$tmpdir/face.mnc",'-invert_transformation');
    $model_face="$tmpdir/face.mnc";
  }

  # create a defacing grid in stx space
  do_cmd('make_random_grid.pl', '--clobber',
         '--mask', $model_face, $model_face, 
         "$tmpdir/grid.mnc",'--amplitude', 
          $amp, '--fwhm', $fwhm);

  do_cmd('itk_morph','--exp','D[2]',$brain_mask,"$tmpdir/brain.mnc",'--clobber');
  do_cmd('itk_morph','--exp','D[1]',$model_face,"$tmpdir/face.mnc" ,'--clobber');

  do_cmd('mincresample',"$tmpdir/brain.mnc","$tmpdir/brain2.mnc",'-like',"$tmpdir/face.mnc",'-nearest','-clobber');
  do_cmd('minccalc','-expression','A[0]==1&&A[1]==0?1:0', "$tmpdir/face.mnc", "$tmpdir/brain2.mnc", "$tmpdir/face2.mnc",'-clobber');
  do_cmd('rm','-f',"$tmpdir/brain.mnc","$tmpdir/brain2.mnc");

  do_cmd('mincresample',"$tmpdir/face2.mnc","$tmpdir/face.mnc",'-like',"$tmpdir/grid.mnc",'-nearest','-clobber');

  do_cmd('mincconcat', '-clobber', '-concat_dimension', 'vector_dimension', '-coordlist',"0,1,2", "$tmpdir/face.mnc","$tmpdir/face.mnc","$tmpdir/face.mnc", "$tmpdir/face2.mnc",'-clobber');
  do_cmd('minccalc', '-expression', 'A[0]*A[1]',"$tmpdir/face2.mnc", "$tmpdir/grid.mnc", "$tmpdir/face_grid.mnc",'-clobber');
  do_cmd('cp',"$tmpdir/face_grid.mnc",$out_grid);
}

deface_volume("$tmpdir/face_grid.mnc",$t1w_xfm,$t1w,"$tmpdir/deface_t1w.mnc");
deface_volume("$tmpdir/face_grid.mnc",$t2w_xfm,$t2w,"$tmpdir/deface_t2w.mnc") if $t2w;
deface_volume("$tmpdir/face_grid.mnc",$pdw_xfm,$pdw,"$tmpdir/deface_pdw.mnc") if $pdw;

$ENV{MINC_COMPRESS}=$compress if $compress;

if($watermark)
{
  watermark("$tmpdir/deface_t1w.mnc",$out_t1w);
  watermark("$tmpdir/deface_t2w.mnc",$out_t2w) if $t2w;
  watermark("$tmpdir/deface_pdw.mnc",$out_pdw) if $pdw;
} else {
  do_cmd('mincreshape',"$tmpdir/deface_t1w.mnc",$out_t1w,'-clobber');
  do_cmd('mincreshape',"$tmpdir/deface_t2w.mnc",$out_t2w,'-clobber') if $t2w;
  do_cmd('mincreshape',"$tmpdir/deface_pdw.mnc",$out_pdw,'-clobber') if $pdw;
}


sub deface_volume {

  my ($deface_grid,$xfm,$scan,$out)=@_;

  do_cmd('xfminvert',$xfm,"$tmpdir/native.xfm",'-clobber');
  do_cmd('uniformize_minc.pl',$scan,"$tmpdir/scan.mnc",'--step',2,'--resample','nearest','--clobber');
  do_cmd('resample_grid',$deface_grid,"$tmpdir/native.xfm","$tmpdir/deface_grid_0.mnc",'--like',"$tmpdir/scan.mnc",'--clobber');

  open XFM,">$tmpdir/deface.xfm" or die;
  print XFM "MNI Transform File\nTransform_Type = Grid_Transform;\nDisplacement_Volume = deface_grid_0.mnc;\n";
  close XFM;

  #my ($lo,$hi)=split(/\n/,`mincstats -q -min -max $scan`);
  my @arg=('mincresample','-nearest','-transform',"$tmpdir/deface.xfm",$scan,$out,'-clobber','-use_input_sampling');#,'-range',$lo,$hi
  push @arg,'-keep_real_range' if $keep_real_range;
  do_cmd(@arg);
  do_cmd('rm','-f',"$tmpdir/native.xfm","$tmpdir/native_grid_0.mnc","$tmpdir/scan.mnc","$tmpdir/deface_grid_0.mnc","$tmpdir/deface.xfm","$tmpdir/deface_grid_0.mnc");

}

sub watermark {

  my ($infile,$outfile)=@_;

  my $posterior = "$model_dir/watermark_posterior.mnc";
  my $inferior  = "$model_dir/watermark_inferior.mnc";
  my $left      = "$model_dir/watermark_left.mnc";

  my %info=minc_info($infile);

  do_cmd('mincreshape','-dimorder','xspace,zspace,yspace',$infile,"$tmpdir/sample.mnc",'-clob');
  reshape_like($posterior,"$tmpdir/sample.mnc","$tmpdir/posterior.mnc");
  reshape_like($inferior ,"$tmpdir/sample.mnc","$tmpdir/inferior.mnc");
  reshape_like($left     ,"$tmpdir/sample.mnc","$tmpdir/left.mnc");

  my $avg=`mincstats -biModalT -q $infile`;
  chomp($avg);
  $avg=int(0.8*$avg);

  do_cmd('mincmath','-nocheck_dimensions','-max',"$tmpdir/posterior.mnc","$tmpdir/inferior.mnc","$tmpdir/left.mnc","$tmpdir/tmp.mnc",'-clobber');
  do_cmd('minccalc','-expression',"A[0]*$avg","$tmpdir/tmp.mnc","$tmpdir/combined.mnc",'-short');

  do_cmd('mincmath','-clobber','-add',"$tmpdir/sample.mnc","$tmpdir/combined.mnc","$tmpdir/out.mnc",'-copy_header','-short','-nocheck_dimensions');

  my @arg=('mincreshape','-dimorder',$info{dimnames},'-clobber',"$tmpdir/out.mnc",$outfile);
  push @arg,'-valid_range',0,4095 unless $no_int_norm;
  do_cmd(@arg);
  
  do_cmd('rm','-f',"$tmpdir/sample.mnc","$tmpdir/combined.mnc","$tmpdir/out.mnc","$tmpdir/posterior.mnc","$tmpdir/inferior.mnc","$tmpdir/left.mnc");
}

sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}
sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}

sub minc_info {
   my ($input)=@_;
   my %info = (
   'dimnames' => undef,
   'xspace'   => undef,
   'yspace'   => undef,
   'zspace'   => undef,
   'xstart' => undef,
   'ystart' => undef,
   'zstart' => undef,
   'xstep' => undef,
   'ystep' => undef,
   'zstep' => undef,
   );   
   ($info{dimnames},
    $info{xspace},$info{yspace},$info{zspace},
    $info{xstart},$info{ystart},$info{zstart},
    $info{xstep},$info{ystep},$info{zstep})= 
    split(/\n/, `mincinfo -vardims image -dimlength xspace -dimlength yspace -dimlength zspace -attvalue xspace:start -attvalue yspace:start -attvalue zspace:start -attvalue xspace:step -attvalue yspace:step -attvalue zspace:step $input`);
    for (values %info) 
    {  
      
      if( /space/ ) 
      { 
        s/\s/,/g; 
      } else  {
        $_*=1.0;
      }
    } #convert to floats
  chop($info{dimnames}); #remove last comma
  #print join(' ',%info);
  return %info;
}

# makes a minc file with the same dimension order and step sign
sub reshape_like {
  my ($in,$sample,$out)=@_;
  my %info=minc_info($sample);
  
  do_cmd('mincreshape',$in,"$tmpdir/tmp.mnc",
    '-dimrange',"xspace=0,$info{xspace}",
    '-dimrange',"yspace=0,$info{yspace}",
    '-dimrange',"zspace=0,$info{zspace}",
    '-dimorder',$info{dimnames},
    '-clobber'
  );
         
  do_cmd('mincreshape',"$tmpdir/tmp.mnc",$out,#"$tmpdir/tmp.mnc",
    '-dimorder',$info{dimnames},
    '-dimsize', 'xspace=-1',
    '-dimsize', 'yspace=-1',
    '-dimsize', 'zspace=-1', 
    $info{xstep}>0?'+xdirection':'-xdirection',
    $info{ystep}>0?'+ydirection':'-ydirection',
    $info{zstep}>0?'+zdirection':'-zdirection',
    '-clobber');
}


#resample like another minc file
sub resample_like {
  my ($in,$sample,$out)=@_;

  my %info=minc_info($sample);
  
  if($info{xstep}<0)
  {
    $info{xstart}+=$info{xstep}*$info{xspace};
    $info{xstep}= -$info{xstep};
  }

  if($info{ystep}<0)
  {
    $info{ystart}+=$info{ystep}*$info{yspace};
    $info{ystep}= -$info{ystep};
  }

  if($info{zstep}<0)
  {
    $info{zstart}+=$info{zstep}*$info{zspace};
    $info{zstep}= -$info{zstep};
  }

  $info{xlen}=$info{xstep}*$info{xspace};
  $info{ylen}=$info{ystep}*$info{yspace};
  $info{zlen}=$info{zstep}*$info{zspace};
  
  my @att=split(/\n/,`mincinfo -attvalue xspace:direction_cosines -attvalue yspace:direction_cosines -attvalue zspace:direction_cosines $sample`);
  
  my @cosx=split(/\s/,$att[0]);
  my @cosy=split(/\s/,$att[1]);
  my @cosz=split(/\s/,$att[2]);
  
  do_cmd('mincreshape',$in,"$tmpdir/tmp.mnc",'-clobber');
  do_cmd('minc_modify_header',
         '-dappend',"xspace:direction_cosines=$cosx[0]",
         '-dappend',"xspace:direction_cosines=$cosx[1]",
         '-dappend',"xspace:direction_cosines=$cosx[2]",

         '-dappend',"yspace:direction_cosines=$cosy[0]",
         '-dappend',"yspace:direction_cosines=$cosy[1]",
         '-dappend',"yspace:direction_cosines=$cosy[2]",

         '-dappend',"zspace:direction_cosines=$cosz[0]",
         '-dappend',"zspace:direction_cosines=$cosz[1]",
         '-dappend',"zspace:direction_cosines=$cosz[2]",

         '-dinsert',"xspace:start=$info{xstart}",
         '-dinsert',"yspace:start=$info{ystart}",
         '-dinsert',"zspace:start=$info{zstart}",
         
         '-dinsert',"xspace:step=$info{xstep}",
         '-dinsert',"yspace:step=$info{ystep}",
         '-dinsert',"zspace:step=$info{zstep}",
         "$tmpdir/tmp.mnc");
  do_cmd('mincresample','-like',$sample,'-clobber','-nearest',"$tmpdir/tmp.mnc",$out);
}

sub correct {
  my ($in,$model,$out)=@_;
  do_cmd("nu_correct", "-clobber", "-iter", 100, "-stop", 0.0001, "-fwhm", 0.1,$in, "$tmpdir/nuc.mnc",'-clobber');
  do_cmd('volume_pol',"$tmpdir/nuc.mnc",$model,'--order',1,'--expfile',"$tmpdir/pol.exp",'--clobber');
  do_cmd('minccalc','-expfile',"$tmpdir/pol.exp","$tmpdir/nuc.mnc",$out,'-clobber');
}

sub fix_sampling {
  my $in_minc=$_[0];
  my $need_fixing=0;
  my $spc;
  foreach $spc(split(/\n/,`mincinfo -attvalue xspace:spacing -attvalue yspace:spacing -attvalue zspace:spacing  $in_minc`))
  {
    $need_fixing=1 if $spc=~/irregular/;
  }
  my $out=$in_minc;

  if($need_fixing)
  {
    $out=$tmpdir.'/'.basename($in_minc,'.gz');
    if($in_minc=~/.gz$/)
    {
      do_cmd("gunzip -c $in_minc >$out");
    } else {
      do_cmd('cp',$in_minc,$out);
    }
    do_cmd('minc_modify_header','-sinsert','xspace:spacing=regular__','-sinsert','zspace:spacing=regular__','-sinsert','yspace:spacing=regular__',$out)
  }
  return $out;
}
