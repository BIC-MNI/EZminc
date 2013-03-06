#!/usr/bin/env perl

##############################################################################
# 
# pipeline_qc.pl
#
# Input:
#
# Output:
#      o a jpeg file with the mask overlaid on the T1 image.
#      o a jpeg file with the T1 image overlaid on the T2 image.
#
#
# Larry Baer, May, 2005
# McConnell Brain Imaging Centre, 
# Montreal Neurological Institute, 
# McGill University
##############################################################################

use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;

my $me = basename($0);
my $verbose = 0;
my $fake = 0;
my $clobber = 0;
my $title;
my $mask;
my $spectral;
my @image_range;
my $spectral_mask;
my @mask_range;
my $cyanred;
my $cyanred_mask;
my $lut;
my $gray_mask;
my $discrete;
my $big;
my $labels_lut;
my $labels;
my $labels_mask;
my $red;
my $green_mask;
my $clamp;

GetOptions(
	   'verbose' => \$verbose,
	   'fake'    => \$fake,
	   'clobber' => \$clobber,
	   'title=s'          => \$title,
	   'mask=s' 	      => \$mask,
	   'spectral'         => \$spectral,
	   'spectral-mask'    => \$spectral_mask,
     'gray-mask'    => \$gray_mask,
	   'image-range=f{2}' => \@image_range,
	   'mask-range=f{2}' => \@mask_range,
     'cyanred'         => \$cyanred,
     'cyanred-mask'    => \$cyanred_mask,
     'lut=s'           => \$lut,
     'discrete'        => \$discrete,
     'big'             => \$big,
     'labels'          => \$labels,
     'labels-mask'     => \$labels_mask,
     'red'             => \$red,
     'green-mask'      => \$green_mask,
     'clamp'           => \$clamp,
	   );

die "Usage: $me  <input.mnc> <output jpeg file> 
[ 
  --title <title> 
  --clobber 
  --verbose 
  --mask <file> 
  --spectral 
  --spectral-mask
  --image-range a b 
  --spectral-mask 
  --gray-mask
  --mask-range a b
  --cyanred 
  --cyanred-mask
  --lut <lut>
  --labels
  --labels-mask
  --dicrete
  --red
  --green-mask
  --clamp
]\n"  if $#ARGV < 1;

#########################
# Get the arguments.
my $infile;

my $out=pop @ARGV;
my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

labels_map();
delete $ENV{MINC_COMPRESS} if $ENV{MINC_COMPRESS};

#print $#image_range,"\n";

foreach $infile(@ARGV)
{
  # make tmpdir

  my $outfile;
  unless($#ARGV) #we should have directory
  {
    $outfile=$out;
  } else {
    $outfile=$out.'/'.basename($infile).'.jpg';
    $title=basename($infile);
  }
  if($clobber || check_file($outfile) )
  {
    my @args=('minclookup','-clobb',#'-range',10,100,
      $infile, "$tmpdir/grey.mnc",'-clobber');

    if($spectral)
    {
      push @args,'-spectral';
    } elsif($cyanred) { 
      push @args,'-lut_string',"0.000 0.8 1.0 1.0;0.125 0.4 0.9 1.0;0.250 0.0 0.6 1.0;\
0.375 0.0 0.2 0.5;0.500 0.0 0.0 0.0;0.625 0.5 0.0 0.0;0.750 1.0 0.4 0.0;\
0.825 1.0 0.8 0.4;1.000 1.0 0.8 0.8";
    } elsif($lut) {
      push @args,'-lookup_table',$lut;
    } elsif($labels) {
      push @args,'-discrete','-lut_string',$labels_lut;
    } elsif($red) {
      push @args,'-lut_string','0.0 0.0 0.0 0.0;1.0 1.0 0.0 0.0';
    } else {
      push @args,'-gray';
    }
  
    push @args,'-discrete' if($discrete);

    if($#image_range>0 && $image_range[0]!=$image_range[1])
    {
      push @args,'-range',$image_range[0],$image_range[1];
    }
    
    do_cmd(@args);
    
    #red mask
    if($mask)
    {
      my @range=split(/\n/,`mincstats -q -min -max $mask`);
      my @args=('minclookup','-clobb',$mask, "$tmpdir/mask.mnc",'-clobber');

      if($spectral_mask)
      {
        push @args,'-spectral';
      } elsif($cyanred_mask) { 
        push @args,'-lut_string',"0.000 0.8 1.0 1.0;0.125 0.4 0.9 1.0;0.250 0.0 0.6 1.0;\
0.375 0.0 0.2 0.5;0.500 0.0 0.0 0.0;0.625 0.5 0.0 0.0;0.750 1.0 0.4 0.0;\
0.825 1.0 0.8 0.4;1.000 1.0 0.8 0.8";
      } elsif($gray_mask) {
        push @args,'-gray';
      } elsif($green_mask) {
        push @args,'-lut_string','0.0 0.0 0.0 0.0;1.0 0.0 1.0 0.0';
      } elsif($labels_mask) {
        push @args,'-discrete','-lut_string',$labels_lut;
      } else {
        push @args,'-lut_string','0.0 0.0 0.0 0.0;1.0 1.0 0.0 0.0';
      }
			
      if($#mask_range>0 && $mask_range[0]!=$mask_range[1])
      {
        push @args,'-range',$mask_range[0],$mask_range[1];
      } else {
				push @args,'-range',$range[0],$range[1];
			}

      do_cmd(@args);
      
      if($clamp)
      {
        do_cmd('minccalc','-byte','-express','clamp(A[0]*1.0+A[1]*0.4,0,1.0)',"$tmpdir/grey.mnc","$tmpdir/mask.mnc","$tmpdir/masked.mnc");
      }else {
        do_cmd('mincmath','-max',"$tmpdir/grey.mnc","$tmpdir/mask.mnc","$tmpdir/masked.mnc");
      }
      do_cmd('mv',"$tmpdir/masked.mnc","$tmpdir/grey.mnc");
    }
    if($big) {
      make_multipane_big("$tmpdir/grey.mnc", $title, $outfile);
    } else {
      make_multipane("$tmpdir/grey.mnc", $title, $outfile);
    }
  }
}


#####################################################################
# sub-routine to make a multipane view 
sub make_multipane{
   my($mncfile, $text, $imgfile, @ext_args) = @_;
   my $smalltilesize = 150;
   my(@args, @mont_args);
   
   # try a .gz if file missing
   $mncfile .= '.gz' if (!-e $mncfile);
   my($xspace,$yspace,$zspace)=split(/\n/,`mincinfo -dimlength xspace -dimlength yspace  -dimlength zspace $mncfile`);
   $xspace=$xspace*1.0;
   $yspace=$yspace*1.0;
   $zspace=$zspace*1.0;
   my $i;
   #z 
   foreach  $i(30.0,35.0,40.0,45.0,50.0,145.0) {
       @args = ('mincpik', '-scale','1','-transverse','-slice',int($zspace*$i/181.0),'-clobber'); #linux:add -clobber
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/T${i}.miff" );
			 do_cmd(@args);
       
       push(@mont_args, "$tmpdir/T${i}.miff");
       }
 #x
  foreach  $i(50,60,70,130,120,110) {
       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$i/181.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$i.miff");
       do_cmd(@args);
       
       push(@mont_args, "$tmpdir/S$i.miff");
       }
 #y
  foreach  $i(60,80,110,120,140,160) {
     @args = ('mincpik','-scale','1','-coronal','-slice',int($yspace*$i/217.0),'-clobber');
     push(@args, @ext_args) if @ext_args;
     push(@args, $mncfile, "$tmpdir/C$i.miff");
     &do_cmd(@args);
     
     push(@mont_args, "$tmpdir/C$i.miff");
     }
 
    # do the montage
  &do_cmd('montage',
          '-tile', '3x6',
          '-background', 'grey10',
          '-geometry', $smalltilesize . 'x' . $smalltilesize . '+1+1',
          @mont_args,
          "$tmpdir/mont.miff");
             
    # Add the title
  @args = ('convert', '-box', 'white');
     #'-font', '7x13', 
     #'-fill', 'white',
     #'-draw', "text 2,15 \"$text\"");
   push(@args,'-draw', "text 2,15 \"$text\"") if $text;
   #push(@args, @more_args) if @more_args;
  &do_cmd(@args,"$tmpdir/mont.miff", $imgfile);
    
 }

#####################################################################
# sub-routine to make a multipane view 
sub make_multipane_big{
   my($mncfile, $text, $imgfile, @ext_args) = @_;
   my $smalltilesize = 170;
   my(@args, @mont_args);

   my $steps=20;
   my $columns=$steps/2;
   my $rows=6;

   # try a .gz if file missing
   $mncfile .= '.gz' if (!-e $mncfile);
   my($xspace,$yspace,$zspace)=split(/\n/,`mincinfo -dimlength xspace -dimlength yspace  -dimlength zspace $mncfile`);
   $xspace=$xspace*1.0;
   $yspace=$yspace*1.0;
   $zspace=$zspace*1.0;
   my $i;
   my $j;
   #z 
   for($j=0;$j<$steps;$j++) {

       $i=int(10+(150.0-10.0)*$j/($steps-1));

       @args = ('mincpik', '-scale','1','-transverse','-slice',int($zspace*$i/181.0),'-clobber'); #linux:add -clobber
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/T${i}.miff" );
			 do_cmd(@args);
       
       push(@mont_args, "$tmpdir/T${i}.miff");
       }
 #x
   for($j=0;$j<($steps/2);$j++) {

       $i=int(28.0+(166.0-28.0)*$j/($steps-1));

       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$i/193.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$i.miff");
       do_cmd(@args);
       
       push(@mont_args, "$tmpdir/S$i.miff");
       }

   for($j=($steps-1);$j>=($steps/2);$j--) {

       $i=int(28.0+(166.0-28.0)*$j/($steps-1));

       @args = ('mincpik','-scale','1','-sagittal','-slice',int($xspace*$i/193.0),'-clobber');
       push(@args, @ext_args) if @ext_args;
       push(@args, $mncfile, "$tmpdir/S$i.miff");
       do_cmd(@args);
       
       push(@mont_args, "$tmpdir/S$i.miff");
       }
 #y
   for($j=0;$j<$steps;$j++) {

     $i=int(25+(195.0-25.0)*$j/($steps-1));
     @args = ('mincpik','-scale','1','-coronal','-slice',int($yspace*$i/217.0),'-clobber');
     push(@args, @ext_args) if @ext_args;
     push(@args, $mncfile, "$tmpdir/C$i.miff");
     &do_cmd(@args);
     
     push(@mont_args, "$tmpdir/C$i.miff");
     }
 
    # do the montage
  &do_cmd('montage',
          '-tile', "${columns}x${rows}",
          '-background', 'grey10',
          '-geometry', $smalltilesize . 'x' . $smalltilesize . '+1+1',
          @mont_args,
          "$tmpdir/mont.miff");
             
    # Add the title
  @args = ('convert', '-box', 'white');
     #'-font', '7x13', 
     #'-fill', 'white',
     #'-draw', "text 2,15 \"$text\"");
   push(@args,'-draw', "text 2,15 \"$text\"") if $text;
   #push(@args, @more_args) if @more_args;
  &do_cmd(@args,"$tmpdir/mont.miff", $imgfile);
    
 }



#####################################################################
sub do_cmd { 
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
      system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  if(-e $_[0])
  {
    warn("${_[0]} exists!\n");
    return 0;
  }    
  return 1;
}

sub labels_map {
 $labels_lut="\
 0 0 0 0;\
 1 1 0 0;\
 2 0 1 0;\
 3 0 0 1;\
 4 0 1 1;\
 5 1 0 1;\
 6 1 1 0;\
 7 0.541176 0.168627 0.886275;\
 8 1 0.0784314 0.576471;\
 9 0.678431 1 0.184314;\
 10 0.12549 0.698039 0.666667;\
 11 0.282353 0.819608 0.8;\
 12 0.627451 0.12549 0.941176;\
 13 1 1 1;\
 14 0.4 0 0;\
 15 0.4 0.2 0;\
 16 0.4 0.4 0;\
 17 0.6 0.4 0;\
 18 0 0.4 0;\
 19 0 0.4 0.2;\
 20 0 0.4 0.4;\
 21 0 0.2 0.4;\
 22 0 0 0.4;\
 23 0.2 0 0.4;\
 24 0.4 0 0.4;\
 25 0.4 0 0.2;\
 26 0.4 0 0;\
 27 0.4 0.2 0;\
 28 0.4 0.4 0;\
 29 0.2 0.4 0;\
 30 0 0.4 0;\
 31 0 0.4 0.2;\
 32 0 0.4 0.4;\
 33 0 0.2 0.4;\
 34 0 0 0.4;\
 35 0.2 0 0.4;\
 36 0.4 0 0.4;\
 37 0.4 0 0.2;\
 38 0.4 0 0;\
 39 0.4 0.2 0;\
 40 0.4 0.4 0;\
 41 0.2 0.4 0;\
 42 0 0.4 0;\
 43 0 0.4 0.2;\
 44 0 0.4 0.4;\
 45 0 0.2 0.4;\
 46 0 0 0.4;\
 47 0.2 0 0.4;\
 48 0.4 0 0.4;\
 49 0.4 0 0.2;\
 50 0.701961 0 0;\
 51 0.701961 0.34902 0;\
 52 0.701961 0.701961 0;\
 53 0.34902 0.701961 0;\
 54 0 0.701961 0;\
 55 0 0.701961 0.34902;\
 56 0 0.701961 0.701961;\
 57 0 0.34902 0.701961;\
 58 0 0 0.701961;\
 59 0.34902 0 0.701961;\
 60 0.701961 0 0.701961;\
 61 0.701961 0 0.34902;\
 62 0.4 0 0;\
 63 0.4 0.2 0;\
 64 0.4 0.4 0;\
 65 0.2 0.4 0;\
 66 0 0.4 0;\
 67 0 0.4 0.2;\
 68 0 0.4 0.4;\
 69 0 0.2 0.4;\
 70 0 0 0.4;\
 71 0.2 0 0.4;\
 72 0.4 0 0.4;\
 73 0.4 0 0.2;\
 74 0.54902 0 0;\
 75 0.54902 0.27451 0;\
 76 0.54902 0.54902 0;\
 77 0.27451 0.54902 0;\
 78 0 0.54902 0;\
 79 0 0.54902 0.27451;\
 80 0 0.54902 0.54902;\
 81 0 0.27451 0.54902;\
 82 0 0 0.54902;\
 83 0.27451 0 0.54902;\
 84 0.54902 0 0.54902;\
 85 0.54902 0 0.27451;\
 86 0.701961 0 0;\
 87 0.701961 0.34902 0;\
 88 0.701961 0.701961 0;\
 89 0.34902 0.701961 0;\
 90 0 0.701961 0;\
 91 0 0.701961 0.34902;\
 92 0 0.701961 0.701961;\
 93 0 0.34902 0.701961;\
 94 0 0 0.701961;\
 95 0.34902 0 0.701961;\
 96 0.701961 0 0.701961;\
 97 0.701961 0 0.34902;\
 98 0.85098 0 0;\
 99 0.85098 0.423529 0;\
 100 0.85098 0.85098 0;\
 101 0.423529 0.85098 0;\
 102 0 0.85098 0;\
 103 0 0.85098 0.423529;\
 104 0 0.85098 0.85098;\
 105 0 0.423529 0.85098;\
 106 0 0 0.85098;\
 107 0.423529 0 0.85098;\
 108 0.85098 0 0.85098;\
 109 0.85098 0 0.423529;\
 110 0.4 0 0;\
 111 0.4 0.2 0;\
 112 0.4 0.4 0;\
 113 0.2 0.4 0;\
 114 0 0.4 0;\
 115 0 0.4 0.2;\
 116 0 0.4 0.4;\
 117 0 0.2 0.4;\
 118 0 0 0.4;\
 119 0.2 0 0.4;\
 120 0.4 0 0.4;\
 121 0.4 0 0.2;\
 122 0.47451 0 0;\
 123 0.47451 0.239216 0;\
 124 0.47451 0.47451 0;\
 125 0.239216 0.47451 0;\
 126 0 0.47451 0;\
 127 0 0.47451 0.239216;\
 128 0 0.47451 0.47451;\
 129 0 0.239216 0.47451;\
 130 0 0 0.47451;\
 131 0.239216 0 0.47451;\
 132 0.47451 0 0.47451;\
 133 0.47451 0 0.239216;\
 134 0.54902 0 0;\
 135 0.54902 0.27451 0;\
 136 0.54902 0.54902 0;\
 137 0.27451 0.54902 0;\
 138 0 0.54902 0;\
 139 0 0.54902 0.27451;\
 140 0 0.54902 0.54902;\
 141 0 0.27451 0.54902;\
 142 0 0 0.54902;\
 143 0.27451 0 0.54902;\
 144 0.54902 0 0.54902;\
 145 0.54902 0 0.27451;\
 146 0.623529 0 0;\
 147 0.623529 0.313725 0;\
 148 0.623529 0.623529 0;\
 149 0.313725 0.623529 0;\
 150 0 0.623529 0;\
 151 0 0.623529 0.313725;\
 152 0 0.623529 0.623529;\
 153 0 0.313725 0.623529;\
 154 0 0 0.623529;\
 155 0.313725 0 0.623529;\
 156 0.623529 0 0.623529;\
 157 0.623529 0 0.313725;\
 158 0.701961 0 0;\
 159 0.701961 0.34902 0;\
 160 0.701961 0.701961 0;\
 161 0.34902 0.701961 0;\
 162 0 0.701961 0;\
 163 0 0.701961 0.34902;\
 164 0 0.701961 0.701961;\
 165 0 0.34902 0.701961;\
 166 0 0 0.701961;\
 167 0.34902 0 0.701961;\
 168 0.701961 0 0.701961;\
 169 0.701961 0 0.34902;\
 170 0.776471 0 0;\
 171 0.776471 0.388235 0;\
 172 0.776471 0.776471 0;\
 173 0.388235 0.776471 0;\
 174 0 0.776471 0;\
 175 0 0.776471 0.388235;\
 176 0 0.776471 0.776471;\
 177 0 0.388235 0.776471;\
 178 0 0 0.776471;\
 179 0.388235 0 0.776471;\
 180 0.776471 0 0.776471;\
 181 0.776471 0 0.388235;\
 182 0.85098 0 0;\
 183 0.85098 0.423529 0;\
 184 0.85098 0.85098 0;\
 185 0.423529 0.85098 0;\
 186 0 0.85098 0;\
 187 0 0.85098 0.423529;\
 188 0 0.85098 0.85098;\
 189 0 0.423529 0.85098;\
 190 0 0 0.85098;\
 191 0.423529 0 0.85098;\
 192 0.85098 0 0.85098;\
 193 0.85098 0 0.423529;\
 194 0.92549 0 0;\
 195 0.92549 0.462745 0;\
 196 0.92549 0.92549 0;\
 197 0.462745 0.92549 0;\
 198 0 0.92549 0;\
 199 0 0.92549 0.462745;\
 200 0 0.92549 0.92549;\
 201 0 0.462745 0.92549;\
 202 0 0 0.92549;\
 203 0.462745 0 0.92549;\
 204 0.92549 0 0.92549;\
 205 0.92549 0 0.462745;\
 206 0.4 0 0;\
 207 0.4 0.2 0;\
 208 0.4 0.4 0;\
 209 0.2 0.4 0;\
 210 0.8 0.4 0.1;\
 211 0 0.4 0.2;\
 212 0 0.4 0.4;\
 213 0 0.2 0.4;\
 214 0 0 0.4;\
 215 0.2 0 0.4;\
 216 0.4 0 0.4;\
 217 0.4 0 0.2;\
 218 0.439216 0 0;\
 219 0.439216 0.219608 0;\
 220 0.439216 0.439216 0;\
 221 0.219608 0.439216 0;\
 222 0 0.439216 0;\
 223 0 0.439216 0.219608;\
 224 0 0.439216 0.439216;\
 225 0 0.219608 0.439216;\
 226 0 0 0.439216;\
 227 0.219608 0 0.439216;\
 228 0.439216 0 0.439216;\
 229 0.439216 0 0.219608;\
 230 0.47451 0 0;\
 231 0.47451 0.239216 0;\
 232 0.47451 0.47451 0;\
 233 0.239216 0.47451 0;\
 234 0 0.47451 0;\
 235 0 0.47451 0.239216;\
 236 0 0.47451 0.47451;\
 237 0 0.239216 0.47451;\
 238 0 0 0.47451;\
 239 0.239216 0 0.47451;\
 240 0.47451 0 0.47451;\
 241 0.47451 0 0.239216;\
 242 0.513725 0 0;\
 243 0.513725 0.254902 0;\
 244 0.513725 0.513725 0;\
 245 0.254902 0.513725 0;\
 246 0 0.513725 0;\
 247 0 0.513725 0.254902;\
 248 0 0.513725 0.513725;\
 249 0 0.254902 0.513725;\
 250 0 0 0.513725;\
 251 0.254902 0 0.513725;\
 252 0.513725 0 0.513725;\
 253 0.513725 0 0.254902;\
 254 0.54902 0 0;\
 255 0 0 0";
}
