#! /usr/bin/env perl
 
# ------------------------------ MNI Header ----------------------------------
#@NAME       : lobe_segment 
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: computes the stereotaxic segmentation/labelling of an 
#              arbitrary MRI volume using its G/W/CSF classification and
#              it non-linear deformation to a pre-labelled model.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Wed Feb 19, 1997, Louis Collins
#@MODIFIED   : Sun Nov 6, 2005, Louis Collins
#@VERSION    : $Id:$
#-----------------------------------------------------------------------------

use warnings "all";
use FindBin;

use Getopt::Tabular;
use MNI::Startup;
use MNI::Spawn;

use MNI::FileUtilities qw(check_output_dirs check_files); # qw(search_directories);
use MNI::PathUtilities qw(split_path); # qw(replace_ext split_path);

# ------------------------------ start here
&Initialize();

if (defined($outputfile) && -e $outputfile) {
    if ($Clobber) {
	unlink($outputfile);
    }
    else {
	    die "$outputfile already exists (use -clobber to overwrite)";
    }
}

                                # prepare the transforms - forward and backward.
                                # and inverses of each.

($def_inv) = &prepare_concatenated_transform($transform,$stx_xfm);

# ------------------------------ Start segmenting!
#
                                # split up the classified volume to get 
                                # grey, white, csf and background voxels.

($grey, $white, $csf, $bkg) = &separate_classes( $classified, $Template );


                                # grey-matter mask

                                # map the standard discrete atlas through
                                # the inverse of the stx deformation and
                                # then use the custom atlas to parcellate
                                # the classified data (identify all grey
                                # struct)
($grey_template)= &build_smallest_template($GreyLimits, $grey, $def_inv);
($first_atlas)  = &customize_atlas($GreyAtlas , $def_inv, $grey_template, "grey_atlas");
($parcel_atlas) = &merge_atlas( $first_atlas,  $grey, "mask", "parcelated_grey");


                                # map the standard lateral ventricles as
                                # well, and then use them to mask the csf
                                # map

($vent )        = &customize_atlas($CSFAtlas ,  $def_inv, $csf,"vent_atlas");
($lat_vent)     = &merge_atlas    ($vent, $csf,         "mask", "parcelated_vent");

                                # now take care of white matter, using the
                                # same method.  I.e. using the grey matter
                                # template grid definition for the target
                                # space, map the standard white matter map
                                # through the inverse of the non-linear
                                # deformation to get a customized white
                                # matter atlas for this subject.

($white_matter)  = &customize_atlas($WhiteAtlas ,  $def_inv, $grey_template,"white_atlas");
($white_atlas)   = &merge_atlas($white_matter, $white, "mask", "parcelated_white");

                                # take special care of the thalamus,
                                # putamen and globus palidus by customizing
                                # the atlas and then masking it with the
                                # CSF

($BG_template)   = &build_template($BG_Limits, $Template, $def_inv);
($customized_BG) = &customize_atlas($Atlas_level50_BG , $def_inv, $BG_template,"level0.5BG");
($customized_BG2)= &merge_atlas    ($customized_BG,    $csf,     "invert mask", "masked_level0.5BG");


                                # now, starting with the parcelated grey
                                # matter atlas, resample it onto the
                                # standard template, then add the lateral
                                # ventricles, the corrected BG and finally
                                # the white matter

($custom_atlas0) = &merge_atlas($parcel_atlas,  $Template,       "like", "grey_like_target");
($custom_atlas1) = &merge_atlas($custom_atlas0, $lat_vent,      "over", "with_vent");
($custom_atlas2) = &merge_atlas($custom_atlas1, $customized_BG2, "over", "with_vent+BG");
($custom_atlas4) = &merge_atlas($custom_atlas2, $white_atlas,    "under","final");


if (defined ($SurfaceObj)) {
   print ("Creating mask...\n");
   ($mask) = build_brain_mask($SurfaceObj, $custom_atlas2);
   print ("done mask.\n");
   ($custom_atlas5) = &merge_atlas($custom_atlas4, $mask, "mask","final_masked");
   &Spawn("cp $custom_atlas5 $outputfile");
} else {
                                # if the cortical surface mask is not used,
                                # then label the brainstem and cerebellum
                                # properly (which means extract the Cerb_l,
                                # Cerb_r and BS from the custom_grey_atlas
                                # and use it to mask the white classified
                                # data.  The result is placed under the
                                # current result).

   ($brain_stem_and_cerebellum_atlas) = &build_cerebellar_and_bs_atlas($first_atlas,"Cerb+BS_atlas");
   ($brain_stem_and_cerebellum)       = &merge_atlas($brain_stem_and_cerebellum_atlas, 
                                                     $white, "mask", "parcelated_Cerb+BS");
   ($custom_atlas5) = &merge_atlas($custom_atlas4, $brain_stem_and_cerebellum,    "under","final+Cerb+BS");
   &Spawn("cp $custom_atlas5 $outputfile");
}

exit;

# ------------------------------ end here!
# ------------------------------
#@NAME       : build_cerebellar_and_bs_atlas
#@INPUT      : name of customized GM atlas
#@OUTPUT     : mask volume cerebellar and brain stem labels only.
#@RETURNS    : the name of the cerebellar and brain stem atlas
#@CREATED    : originally: 1.19.98  Louis
#@MODIFIED   : 


sub build_cerebellar_and_bs_atlas {
   my ($atlas, $str) = @_ ;
   my $filebase = (&split_path($atlas))[1];

   my $mask = "${TmpDir}/${filebase}_${str}.mnc";

   if (-e $mask) {
      print "$mask exists already, skipping this build..\n";
   }
   else {

      &Spawn("$MincLookup $atlas $mask -clob -discrete -lut_string '76 76; 67 67; 20 20'");
   }
   ($mask);
}

# ------------------------------
#@NAME       : build_brain_mask
#@INPUT      : name of extracted cortical surface
#@OUTPUT     : mask volume with 1's inside and 0's outside the surface
#@RETURNS    : the name of the brain mask
#@CREATED    : originally: 5.21.97 Louis
#@MODIFIED   : 


sub build_brain_mask {
   my ($surf, $template) = @_ ;
   my $filebase = (&split_path($surf))[1];

   my $mask = "${TmpDir}/brainmask_${filebase}.mnc";

   if (-e $mask) {
      print "$mask exists already, skipping this build..\n";
   }
   else {
      my $temp1 = "${TmpDir}/brainmask_template1_${filebase}.mnc";
      my $temp2 = "${TmpDir}/brainmask_template2_${filebase}.mnc";

                                # need to build a volume of 1's

      &Spawn("$Mincmath -scale -const2 0.0 1.0 $template $temp1");
      &Spawn("$MincReshape -image_range 0 255 -valid_range 0 255 $temp1 $temp2");
      unlink $temp1;
      &Spawn("$SurfaceMask $temp2 $surf $temp1");
      unlink $temp2;
      &Spawn("$Dilate $temp1 $mask 1 26 1");
      unlink $temp1;

   }

   ($mask);
}

# ------------------------------
#@NAME       : merge_atlas
#@INPUT      : name of base atlas, name of second atlas, name of operation, name of result
#@OUTPUT     : customized atlas resampled on standard template
#@RETURNS    : the name of the customized atlas
#@CREATED    : originally: 2.11.97 Louis
#@MODIFIED   : 


sub merge_atlas {
  my ($first, $second, $operation, $name) = @_ ;

  my $custom = "${TmpDir}/custom_${name}.mnc";

  my $tmp    = "${TmpDir}/custom_${name}_tmp.mnc";

  if (-e $custom) {
    print "file $custom exists already\n";
  }
  else {

    if ($operation eq "like") {
      &Spawn("$AutoCrop  $first -from $second $tmp");      
    }
    elsif ($operation eq "invert mask") {
      my $tmp2_rsl   = "${TmpDir}/custom_rsl2.mnc";
      &Spawn("$MincResample -clob -nearest -like $first $second $tmp2_rsl");
      &Spawn("$MincCalc -expression 'A[1] > 0 ? 0.0 : A[0]' $first $tmp2_rsl $tmp");      
      unlink $tmp2_rsl;
    }
    elsif ($operation eq "mask") {
      my $tmp2_rsl   = "${TmpDir}/custom_rsl2.mnc";
      &Spawn("$MincResample -clob -nearest -like $first $second $tmp2_rsl");
      &Spawn("$MincCalc -expression 'A[1] > 0 ? A[0] : 0.0' $first $tmp2_rsl $tmp");      
      unlink $tmp2_rsl;
    }
    elsif ($operation eq "add") {
      &Spawn("$Mincmath -add $first $second $tmp");      
    }
    elsif ($operation eq "over") {
      my $tmp_masked = &merge_atlas($first, $second, "invert mask", "tmpmasked");
      my $tmp_rsl    = "${TmpDir}/custom_rsl.mnc";
      my $tmp2_rsl   = "${TmpDir}/custom_rsl2.mnc";
      &Spawn("$MincResample -nearest -like $first $second $tmp2_rsl");
      &Spawn("$MincReshape -image_range 0 255 -valid_range 0 255 $tmp2_rsl $tmp_rsl");
      my $tmp = &merge_atlas($tmp_masked, $tmp_rsl, "add", "${name}_tmp");
      unlink ($tmp_masked, $tmp_rsl, $tmp2_rsl);
    }
    elsif ($operation eq "under") {
      my $tmp_masked = &merge_atlas($second, $first, "invert mask", "tmpmasked");
      my $tmp_rsl    = "${TmpDir}/custom_rsl.mnc";
      my $tmp2_rsl    = "${TmpDir}/custom_rsl2.mnc";
      &Spawn("$MincResample -nearest -like $first $tmp_masked $tmp2_rsl");
      &Spawn("$MincReshape -image_range 0 255 -valid_range 0 255 $tmp2_rsl $tmp_rsl");
      my $tmp = &merge_atlas($first, $tmp_rsl, "add", "${name}_tmp");
      unlink ($tmp_masked, $tmp_rsl, $tmp2_rsl);
    }
    else {
       die "Undefined operation <$operation> in sub merge_atlas\n";
    }

    &Spawn("$MincReshape -image_range 0 255 -valid_range 0 255 $tmp $custom");
    unlink $tmp;
  }
  ($custom)
}

# ------------------------------
#@NAME       : build_template
#@INPUT      : vol to be deformed, vol with target template sampling, deformation
#@OUTPUT     : template volume that will be in a mincresample -like 
#@RETURNS    : the name of the customized template
#@CREATED    : originally:5.21.97 Louis
#@MODIFIED   : 


sub build_template {
  my ($limits, $template, $xfm) = @_ ;

  my ($tags, $result, $output);

  my $filebase = (&split_path($limits))[1];
  my $custom = "${TmpDir}/${filebase}_custom_template.mnc";
  if (-e $custom) {
    print "file $custom exists already\n";
  } else {
     $tags = build_minimax_tags($limits);
     &Spawn("$AutoCrop -byte -isoexpand 5% $template $custom -from $tags");
  }
  
  ($custom)
}

# ------------------------------
#@NAME       : build_smallest_template
#@INPUT      : tags defining vol to be deformed, 
#              target template volume sampling
#              deformation xfm to be applied to tag points to 
#                  get into the target space.
#@OUTPUT :     the smallest template volume containing the object defined by
#              the points while masked by the template volume.  the result will 
#              be used in a mincresample -like
#@RETURNS    : the name of the customized template
#@CREATED    : originally: 3.25.99 Louis
#@MODIFIED   : 


sub build_smallest_template {
  my ($limits, $template_mask, $xfm) = @_ ;

  my ($tags, $result, $output);

  my $filebase = (&split_path($limits))[1];
  my $custom = "${TmpDir}/${filebase}_custom_template.mnc";
  if (-e $custom) {
    print "file $custom exists already\n";
  }
  else {
                                # define the tags in the target space.
     $tags = build_minimax_tags($limits);
                                # build a volume mask representing the object
     my $custom1 = "${TmpDir}/${filebase}_custom_template1.mnc";
     my $custom2 = "${TmpDir}/${filebase}_custom_template2.mnc";
     my $custom3 = "${TmpDir}/${filebase}_custom_template3.mnc";
     &Spawn("$AutoCrop -isoexpand 5% $template_mask $custom1 -from $tags");
     &Spawn("$Mincmath -clob -gt -constant -10e20 $custom1 $custom2");

                                # mask this volume with the template
     &Spawn("$Mincmath -clob -gt -constant -10e20 $template_mask $custom1");
     &Spawn("$MincResample -clob $custom2 $custom3 -like $custom1 -fillvalue 0");
     &Spawn("$Mincmath -clob -zero -and $custom1 $custom3 $custom2");

#    &Spawn("$AutoCrop -clob -bbox  $custom2 $custom2 $custom");
     &Spawn("$AutoCrop -clob $custom2 $custom");

     unlink $custom1, $custom2, $custom3;
  }
  
  ($custom)
}

sub build_minimax_tags {
   my ($tags) = @_;
   my $filebase = (&split_path($tags))[1];
   my $custom = "${TmpDir}/${filebase}_minmax.tag";
   
   my ( $minx, $maxx, $miny, $maxy, $minz, $maxz );


   open (TAG, "$Stats_tag_file $tags |");

   while (<TAG>) {
      $line = $_;
      chomp $line;
      if (/^X/) {
         ($minx, $maxx) = (split(" ",$line))[2,6];
      } 
      else {
         if (/^Y/) {
            ($miny, $maxy) = (split(" ",$line))[2,6];
         } 
         else {
            if (/^Z/) {
               ($minz, $maxz) = (split(" ",$line))[2,6];
            } 
         }
      }
   }
   close TAG;


   open (TAG, ">$custom");
   print TAG "MNI Tag Point File\n";
   print TAG "Volumes = 1;\n";
   print TAG "% temp minimax tag file;\n\nPoints =\n";

   print TAG "$minx $miny $minz \"\"\n";
   print TAG "$minx $miny $maxz \"\"\n";
   print TAG "$minx $maxy $minz \"\"\n";
   print TAG "$minx $maxy $maxz \"\"\n";
   print TAG "$maxx $miny $minz \"\"\n";
   print TAG "$maxx $miny $maxz \"\"\n";
   print TAG "$maxx $maxy $minz \"\"\n";
   print TAG "$maxx $maxy $maxz \"\";\n";

   close TAG;

   ($custom);
}


# ------------------------------
#@NAME       : customize_atlas
#@INPUT      : atlas to be customized + deformation field to use
#@OUTPUT     : customized atlas resampled on standard template
#@RETURNS    : the name of the customized atlas
#@CREATED    : originally: 2.11.97 Louis
#@MODIFIED   : 


sub customize_atlas {
  my ($atlas, $xfm, $template, $name) = @_ ;

  my $custom = "${TmpDir}/custom_${name}.mnc";

  if (-e $custom) {
    print "file $custom exists already\n";
  }
  else {
    &Spawn("$MincResample $atlas $custom -like $template -transformation $xfm -nearest -keep_real_range");
  }
  ($custom)
}


# ------------------------------
#@NAME       : segment
#@INPUT      : $tissue volume, $structure mask, label value
#@OUTPUT     : segmented structure with label
#@RETURNS    : the name of the segmented structure
#@DESCRIPTION: 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : originally: 10.21.96 Louis
#@MODIFIED   : 

sub segment  {
    my ($file,$mask,$label) = @_;

    my ($tmp, $tmp2, $labelled_struct);
    
    $labelled_struct = $file;
    $labelled_struct =~ s/\.mnc/\_$label.mnc/;
    $tmp = $labelled_struct;
    $tmp =~ s/\.mnc/\_tmp.mnc/g;
    $tmp2 = $labelled_struct;
    $tmp2 =~ s/\.mnc/\_tmp2.mnc/g;

    print "labelled: $labelled_struct \n";  


    if (-e $labelled_struct) {
      print "file      $labelled_struct already exists\n";
    }
    else {
      if (-e "$tmp") {
	print ("file $tmp already exists\n");
      }
      else {
	&Spawn("$Mincmath $file $mask -mult $tmp");
      }
    
      &Spawn("$Mincmath $tmp $tmp2 -seg -const2 0.1 100000");
      &Spawn("$Mincmath $tmp2 $labelled_struct -mult -constant $label");

      unlink ($tmp,$tmp2);
    }

    ($labelled_struct);
    
}

# ------------------------------
#@NAME       : separate classes
#@INPUT      : the name of the classified minc file
#@OUTPUT     : three separate volume for grey, white and csf volumes
#@RETURNS    : the name of the subsampled output minc file
#@DESCRIPTION: subsample the input file to a 2x2x2 grid
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : originally: 10.21.96 Louis
#@MODIFIED   : 

sub separate_classes  {
    my ($file,$template) = @_;

    my ($grey, $white, $csf );
    my ($llimit, $ulimit);
				# get the file basename, without directory
    $filebase = (&split_path($file))[1];

    $temp  = "${TmpDir}/${filebase}_like_template.mnc";
    $grey  = "${TmpDir}/${filebase}_grey.mnc";
    $white = "${TmpDir}/${filebase}_white.mnc";
    $csf   = "${TmpDir}/${filebase}_csf.mnc";
    $bkg   = "${TmpDir}/${filebase}_bkg.mnc";

    
    if (-e "$temp") {
	print ("file $temp already exists\n");
    } else {
	&Spawn("$MincResample $file $temp -like $template -nearest -clob");
    }

    if (-e "$grey") {
	print ("file $grey already exists\n");
    } else {
        $llimit = $GWCB[0] - 0.5;
        $ulimit = $GWCB[0] + 0.5;
	&Spawn("$Mincmath -byte $temp $grey -seg -const2 $llimit $ulimit");
    }
    
    if (-e "$white") {
	print ("file $white already exists\n");
    } else {
        $llimit = $GWCB[1] - 0.5;
        $ulimit = $GWCB[1] + 0.5;
	&Spawn("$Mincmath -byte $temp $white -seg -const2 $llimit $ulimit");
    }
    
    if (-e "$csf") {
	print ("file $csf already exists\n");
    } else {
        $llimit = $GWCB[2] - 0.5;
        $ulimit = $GWCB[2] + 0.5;
	&Spawn("$Mincmath -byte $temp $csf -seg -const2 $llimit $ulimit");
    }

    if (-e "$bkg") {
	print ("file $bkg already exists\n");
    } else {
        $llimit = $GWCB[3] - 0.5;
        $ulimit = $GWCB[3] + 0.5;
	&Spawn("$Mincmath -byte $temp $bkg -seg -const2 $llimit $ulimit");
    }

    ($grey, $white, $csf, $bkg );
    
}
# ------------------------------
#@NAME       : preprocess
#@INPUT      : the name of the input minc file
#@OUTPUT     : whatever $Autocrop spits out
#@RETURNS    : the name of the subsampled output minc file
#@DESCRIPTION: subsample the input file to a 2x2x2 grid
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : originally: 10.21.96 Louis
#@MODIFIED   : 

sub preprocess {
    my ($file) = @_;

    my ($result);

				# get the file basename, without directory
    $filebase = (&split_path($file))[1];

    $result = "${TmpDir}/${filebase}_iso2.mnc";

				# build the 2x2x2 mm^3 volume
    if (-e "$result") {
	print ("file $result already exists\n");
    }
    else {
	&Spawn("$AutoCrop $file $result -isostep 2");
    }
    
    ($result);
}

sub prepare_transforms {

  my ($xfm) = @_;

  my ($forward, $forward_rev, 
      $inverse, $inverse_rev,
      $def, $def_inv,
      $filebase);

  $filebase = (&split_path($xfm))[1];

  $forward         = $xfm;
  $forward_rev     = "${TmpDir}/${filebase}_rev.xfm";
  $inverse         = "${TmpDir}/${filebase}_inv.xfm";
  $inverse_rev     = "${TmpDir}/${filebase}_inv_rev.xfm";
  $def             = "${TmpDir}/${filebase}_def_only.xfm";
  $def_inv         = "${TmpDir}/${filebase}_def_only_inv.xfm";
  
  if (-e $forward_rev) {
    print "file $forward_rev already exists.\n";
  }
  else {
    &Spawn("reversedef $xfm $forward_rev");
  }

  if (-e $inverse) {
    print "file $inverse already exists.\n";
  }
  else {
    &Spawn("$XfmTool  invert:${xfm} $inverse");
  }

  if (-e $inverse_rev) {
    print "file $inverse_rev already exists.\n";
  }
  else {
    &Spawn("$XfmTool  invert:${forward_rev} $inverse_rev");
  }

  if ( &is_there_a_grid_transform($forward) ) {

    if (-e $def) {
      print "file $def already exists.\n";
    }
    else {
      &Spawn("$XfmTool  \'extract(2):${forward_rev}\' $def");
    }

    if (-e $def_inv) {
      print "file $def_inv already exists.\n";
    }
    else {
      &Spawn("$XfmTool  \'extract(1):${inverse}\' $def_inv");
    }
  }
  else {
    &Spawn("$Param2Xfm $def");      # make two identity xforms
    &Spawn("$Param2Xfm $def_inv");
  }

  ($forward, $forward_rev, $inverse, $inverse_rev, $def, $def_inv); 
}

sub prepare_concatenated_transform {

  my ($xfm,$stx) = @_;

  my ($forward, 
      $inverse, 
      $def, $def_inv,
      $filebase);

  $filebase = (&split_path($xfm))[1];

  $inverse_stx     = "${TmpDir}/${filebase}_inv_stx.xfm";
  
  if (-e $inverse_stx) {
    print "file $inverse_stx already exists.\n";
  }
  else {
    &Spawn("$XfmTool  invert:${xfm} $stx $inverse_stx");
  }

  ($inverse_stx); 
}

sub is_there_a_grid_transform {
  my ($file) = @_;

  my $res = 0;

  open(FOO, "grep Grid_Transform ${file}|");
  while (<FOO>) {
    if (/Grid_Transform/) {
      $res = 1;
    }
  }
  close(FOO);

  ($res);
}




# --------------------------------------------
sub Initialize
{
   $Version = "1.5";
   $LongVersion = "version ${Version}: slightly tested perl code. Beware!";

   &self_announce if $Verbose && ! -t "STDOUT";
   MNI::Spawn::SetOptions (err_action => "fatal");

   my $usage = <<USAGE;
Usage: $ProgramName [options] to_model_nl.xfm to_tal.xfm  class.mnc out_labels.mnc
       $ProgramName -help for details

USAGE

   my $help = <<HELP;
Help text:

$ProgramName will use the input transformation to define the mapping
between the input MRI volume and the atlas labels to segment lobes and
BG stuctures.  This script takes three arguments:

       1- the name of the non-linear native to model transformation
       2- the name of the linear stereotaxic transformation       
       3- the name of the stereotaxic classified data (one volume
          containing grey, white, csf and background.)
       4- the name of the output label volume (in stereotaxic space).

The corresponding MRI volume is segmented in stereotaxic space in a two
step process.  First the inverse of the non-linear to_model transformation
is concatenated with the linear stereotaxic transformation.  The resulting
transformation is used customize a STANDARD ATLAS (described below) for the
particular subject by resampling the ATLAS through inverse transformation
onto the stereotaxic linearly transformed MRI volume.  This essentially
achieves structure identification, however the fine structure borders
(e.g., grey/white interface) are not usually well defined by this
procedure.

In order to refine the segmentation, structures of the customized atlas are
used to parcelate the stereotaxic classifed data. This step allows the grey
matter volume to be split into different lobes, for example.

After masking, some structures are not always properly segmented since they
do not contain a single pure tissue class: the medial third of the thalamus
is usually classified as grey matter while the rest of it is classified as
white.  In order to fix these segmentations, the customized max-proba atlas
for the basal ganglia is masked with the CSF volume and the result is
merged with the previous grey-matter structure labels.

HELP


  $Atlas_source = "icbm152-lobes-v1.1"; # use the standard ICBM atlas by default

  @GWCB = ( 2, 3, 1, 0 );

  @ArgInfo =
    (@DefaultArgs,
     ["User definable flags", "section"],
     ["-version", "call", undef, \&print_version,
      "print version and quit"],
     ["-surface_mask", "string", 1, \$SurfaceObj,
        "specify and use cortical surface mask (.obj)"],
     ["-template", "string", 1, \$Template,
        "define a target template (.mnc)"],
     ["-modeldir", "string", 1, \$Modeldir,
        "define the directory for the model"],
     ["-gwcb", "integer", 4, \@GWCB,
        "list of indices for grey white csf and background (no commas)"]     
    );
   
   &Getopt::Tabular::SetHelp($help, $usage);
   
   my (@argv) = @ARGV;
   &GetOptions(\@ArgInfo, \@argv) || exit;

   if (@argv != 4) {
     die $usage;
   }
  
   # atlas volumes used to drive segmentation
   if (! $Modeldir) {
     $Modeldir         = "$FindBin::Bin/../share/ANIMAL_INSECT/" . $Atlas_source ;
   }

   $GreyAtlas        = "${Modeldir}/AtlasGrey.mnc";
   $WhiteAtlas       = "${Modeldir}/AtlasWhite.mnc";
   $CSFAtlas         = "${Modeldir}/atlas_csf.mnc";
   $Atlas_level50_BG = "${Modeldir}/atlas_level0.5_BG.mnc";
   $GreyLimits       = "${Modeldir}/limits_grey.tag";
   $VenLimits        = "${Modeldir}/limits_ven_manual.tag";
   $BG_Limits        = "${Modeldir}/limits_level0.5_BG.tag";
 
   @model_files = ($GreyAtlas, $WhiteAtlas, $CSFAtlas, $Atlas_level50_BG,
                   $GreyLimits, $VenLimits, $BG_Limits );

   &check_files(@model_files) || die 'Missing model files';

                                # (spatial definition) template files

   if (! $Template) {
      $Template         =  "$FindBin::Bin/../share/ICBM/" . "icbm_template_1.00mm.mnc";
   }
   &check_files( $Template ) || die 'Missing template file';

   # Look for required programs

   $Dilate       = "dilate_volume";
   $SurfaceMask  = "surface_mask2";
   $AutoCrop     = "autocrop";
   $Mincmath     = "mincmath";
   $XfmTool      = "xfmtool";
   $MincResample = "mincresample";
   $MincReshape  = "mincreshape";
   $MincLookup   = "minclookup";
   $MincCalc     = "minccalc";
   $Stats_tag_file = "stats_tag_file";

   my @programs =
      qw/cp dilate_volume surface_mask2 autocrop
         mincmath xfmtool mincresample
         mincreshape minclookup minccalc stats_tag_file/;
   RegisterPrograms(\@programs);

   # They were found, so add options according to
   # $Debug, $Verbose, $Clobber flags
   
   my ($debug, $verbose, $clobber);
   $debug   = ($Debug)   ? " -debug"   : "";
   $verbose = ($Verbose) ? ""          : " -quiet";
   $clobber = ($Clobber) ? " -clobber" : "";
   
   $AutoCrop     .= "$verbose$clobber";
   $MincResample .= "$verbose$clobber";
   $MincReshape  .= "$verbose$clobber";
   $MincCalc     .= "$debug$verbose$clobber";
   $Mincmath     .= "$debug$verbose$clobber";
   
   ($transform, $stx_xfm, $classified, $outputfile) = @argv;

   die "$transform does not exist.\n"  unless (-e $transform || -e "${transform}.gz") ;
   die "$stx_xfm does not exist.\n"    unless (-e $stx_xfm || -e "${stx_xfm}.gz");
   die "$classified does not exist.\n" unless (-e $classified || -e "${classified}.gz");
 

   &check_output_dirs($TmpDir);
}

sub print_version  {
    die "Program $ProgramName, built from:\n$LongVersion\n";
}
