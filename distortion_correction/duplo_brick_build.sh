#! /bin/sh

#dimensions="-start -50 -50 -50 -nelements 100 100 100  -step 1 1 1"
dimensions="-start -40 -25 -25 -nelements 320 200 200  -step 0.25 0.25 0.25"
#dimensions="-start -40 -25 -25 -nelements 750 500 500  -step 0.1 0.1 0.1"

tempdir=$(mktemp -d --tmpdir)
trap "rm -rf $tempdir" 0 1 2 15

calc="-max_buffer_size_in_kb 10000000"

#unset MINC_FORCE_V2
#unset MINC_COMPRESS
make_phantom -rectangle -width 62.75 31  18.5  $dimensions $tempdir/rect_outer.mnc -center 0 0 0.5 -clob
make_phantom -rectangle -width 60.5  28.5 17.5 $dimensions $tempdir/rect_inner.mnc -center 0 0 1.5 -clob
mincmath -sub $tempdir/rect_outer.mnc $tempdir/rect_inner.mnc $tempdir/rect_shell.mnc -clob $calc

indent_top_bottom="-rectangle -width 1 2 18.5"
indent_left_right="-rectangle -width 2 1 18.5"
#make_phantom $indent $dimensions $tempdir/rect_indent.mnc        -center 0    0 0.5
#make_phantom $indent $dimensions $tempdir/rect_indent_8_8.mnc    -center 8    8 0.5
#make_phantom $indent $dimensions $tempdir/rect_indent_14_8.mnc   -center 14   8 0.5

#top row
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_8_14.mnc    -center   8  14 0.5
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_24_14.mnc   -center  24  14 0.5
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_-8_14.mnc   -center  -8  14 0.5
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_-24_14.mnc  -center -24  14 0.5
#bottom row
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_8_-14.mnc   -center   8 -14 0.5
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_24_-14.mnc  -center  24 -14 0.5
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_-8_-14.mnc  -center  -8 -14 0.5
make_phantom $indent_top_bottom $dimensions $tempdir/rect_indent_-24_-14.mnc -center -24 -14 0.5

#right column
make_phantom $indent_left_right $dimensions $tempdir/rect_indent_30_8.mnc    -center   30  8 0.5
make_phantom $indent_left_right $dimensions $tempdir/rect_indent_30_-8.mnc   -center   30 -8 0.5
#left column
make_phantom $indent_left_right $dimensions $tempdir/rect_indent_-30_8.mnc   -center  -30  8 0.5
make_phantom $indent_left_right $dimensions $tempdir/rect_indent_-30_-8.mnc  -center  -30 -8 0.5

mincmath -max $tempdir/rect_indent_*   $tempdir/all_indents.mnc -clobber $calc
mincmath -max $tempdir/all_indents.mnc $tempdir/rect_shell.mnc $tempdir/rect_carcas.mnc $calc


param2xfm -translation 16 0 0  $tempdir/move_16.xfm
param2xfm -translation -16 0 0 $tempdir/move_-16.xfm

make_phantom -ellipse -xwidth 13 -ywidth 13 -zwidth 2000 $dimensions $tempdir/ellipse.mnc -center 0 0 0 -clobber
make_phantom -ellipse -xwidth 11 -ywidth 11 -zwidth 2000 $dimensions $tempdir/ellipse_inner.mnc -center 0 0 0 -clobber
mincmath -sub $tempdir/ellipse.mnc $tempdir/ellipse_inner.mnc $tempdir/cylinder.mnc $calc

make_phantom -rectangle -xwidth 20 -ywidth 20 -zwidth 16 $dimensions $tempdir/ellipse_cut.mnc   -center 0 0 0 -clobber
minccalc -express 'A[1]>0.5?A[0]:0' $tempdir/cylinder.mnc $tempdir/ellipse_cut.mnc $tempdir/center_cylinder.mnc -clob  $calc

mincresample $tempdir/center_cylinder.mnc $tempdir/cyl_16.mnc  -transform  $tempdir/move_16.xfm -nearest -use_input_sampling -clob
mincresample $tempdir/center_cylinder.mnc $tempdir/cyl_-16.mnc -transform $tempdir/move_-16.xfm -nearest -use_input_sampling -clob
mincmath -max $tempdir/center_cylinder.mnc $tempdir/cyl_-16.mnc $tempdir/cyl_16.mnc  $tempdir/cylinders.mnc -clob $calc

# param2xfm -translation 0 0 -2 $tempdir/move_cyl.xfm
# mincresample -nearest -use_input_sampling $tempdir/cylinders.mnc $tempdir/cylinders2.mnc -transform $tempdir/move_cyl.xfm

make_phantom -rectangle -xwidth 0.5 -ywidth 8 -zwidth 12.50 $dimensions $tempdir/wall1.mnc -center 0  10.5 -1.50 -clobber
make_phantom -rectangle -xwidth 0.5 -ywidth 8 -zwidth 12.50 $dimensions $tempdir/wall2.mnc -center 0 -10.5 -1.50 -clobber

mincmath -max $tempdir/rect_carcas.mnc $tempdir/cylinders.mnc $tempdir/wall1.mnc $tempdir/wall2.mnc $tempdir/lower_brick.mnc -clob $calc
make_phantom -ellipse -xwidth 9 -ywidth 9 -zwidth 10000 $dimensions $tempdir/stud.mnc -center 0 0 0 
make_phantom -ellipse -xwidth 7 -ywidth 7 -zwidth 10000 $dimensions $tempdir/stud_inner.mnc -center 0 0 0 
mincmath -sub $tempdir/stud.mnc $tempdir/stud_inner.mnc $tempdir/stud_walls.mnc $calc

#                                    v stud height
make_phantom -rectangle -width 10 10 5.0 $dimensions $tempdir/stud_mask.mnc -center 0 0 -11.5 -clob

minccalc -express 'A[1]>0.5?A[0]:0' $tempdir/stud_walls.mnc $tempdir/stud_mask.mnc $tempdir/stud_c.mnc  $calc


param2xfm -translation 8 8 0  $tempdir/move_8_8.xfm
param2xfm -translation -8 8 0  $tempdir/move_-8_8.xfm
param2xfm -translation -8 -8 0 $tempdir/move_-8_-8.xfm
param2xfm -translation 8 -8 0  $tempdir/move_8_-8.xfm

mincresample -nearest -use_input_sampling $tempdir/stud_c.mnc -transform $tempdir/move_8_8.xfm   $tempdir/stud_8_8.mnc
mincresample -nearest -use_input_sampling $tempdir/stud_c.mnc -transform $tempdir/move_-8_8.xfm  $tempdir/stud_-8_8.mnc
mincresample -nearest -use_input_sampling $tempdir/stud_c.mnc -transform $tempdir/move_-8_-8.xfm $tempdir/stud_-8_-8.mnc
mincresample -nearest -use_input_sampling $tempdir/stud_c.mnc -transform $tempdir/move_8_-8.xfm  $tempdir/stud_8_-8.mnc
mincmath -max $tempdir/stud_-8_-8.mnc $tempdir/stud_-8_8.mnc $tempdir/stud_8_-8.mnc $tempdir/stud_8_8.mnc $tempdir/studs.mnc $calc
mincresample $tempdir/studs.mnc -nearest $tempdir/studs_right.mnc -transform $tempdir/move_16.xfm -use_input_sampling
mincresample $tempdir/studs.mnc -nearest $tempdir/studs_left.mnc -transform  $tempdir/move_-16.xfm -use_input_sampling
mincmath -max $tempdir/studs_right.mnc $tempdir/studs_left.mnc $tempdir/studs_all.mnc $calc
make_phantom -rectangle -xwidth 3 -ywidth 0.5 -zwidth 12 $tempdir/wall_small1.mnc -center 8 0 -2 -clobber $dimensions
make_phantom -rectangle -xwidth 3 -ywidth 0.5 -zwidth 12 $tempdir/wall_small2.mnc -center -8 0 -2 -clobber $dimensions
mincmath -max $tempdir/studs_all.mnc $tempdir/lower_brick.mnc $tempdir/wall_small1.mnc $tempdir/wall_small2.mnc $tempdir/lego_brick.mnc $calc


export MINC_FORCE_V2=1
export MINC_COMPRESS=4


minccalc -expression 'clamp(A[0],0,1)*100' -byte $tempdir/lego_brick.mnc $tempdir/lego_brick_ideal_clamp.mnc      -clob $calc
minccalc -expression 'clamp(A[0],0,1)' -byte $tempdir/rect_outer.mnc $tempdir/lego_brick_ideal_mask_clamp.mnc -clob $calc

# rotate 180 degrees
param2xfm -rotations 180 0 0 $tempdir/rot_x180.xfm

mincresample -use_input -nearest -labels \
 -transform $tempdir/rot_x180.xfm \
 $tempdir/lego_brick_ideal_mask_clamp.mnc \
 duplo_brick_ideal_mask.mnc

mincresample -use_input -nearest -transform  $tempdir/rot_x180.xfm \
    $tempdir/lego_brick_ideal_clamp.mnc \
    duplo_brick_ideal.mnc 

#param2xfm -translation 0 0 18 move_18.xfm 
#mincresample -transform move_18.xfm -use_input_sampling lego_brick_ideal.mnc lego_brick_ideal_moved.mnc -nearest
