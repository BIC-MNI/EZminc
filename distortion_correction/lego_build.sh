#! /bin/sh

#dimensions="-start -50 -50 -50 -nelements 100 100 100  -step 1 1 1"
dimensions="-start -40 -25 -25 -nelements 320 200 200  -step 0.25 0.25 0.25"
calc="-max_buffer_size_in_kb 10000000"
mkdir -p tmp
#unset MINC_FORCE_V2
#unset MINC_COMPRESS
make_phantom -rectangle -width 62.75 31  18.5 $dimensions tmp/rect_outer.mnc -center 0 0 0.5 -clob
make_phantom -rectangle -width 60.5  28.5 17   $dimensions tmp/rect_inner.mnc -center 0 0 1 -clob
mincmath -sub tmp/rect_outer.mnc tmp/rect_inner.mnc tmp/rect_shell.mnc -clob $calc

make_phantom -rectangle -width 2 2 18.5 $dimensions tmp/rect_indent.mnc -center 0 0 0.5
param2xfm -translation 8 8 0 tmp/move_8_8.xfm
mincresample tmp/rect_indent.mnc tmp/rect_indent_8_9.mnc -transform tmp/move_8_8.xfm -nearest -use_input_sampling
rm tmp/rect_indent_8_9.mnc
param2xfm -translation 14 8 0 tmp/move_15_8.xfm
mincresample tmp/rect_indent.mnc tmp/rect_indent_15_8.mnc -transform tmp/move_15_8.xfm -nearest -use_input_sampling
rm tmp/rect_indent_15_8.mnc tmp/move_15_8.xfm
param2xfm -translation 8 14 0 tmp/move_8_15.xfm
mincresample tmp/rect_indent.mnc tmp/rect_indent_8_15.mnc -transform tmp/move_8_15.xfm -nearest -use_input_sampling
param2xfm -translation 24 14 0 tmp/move_24_15.xfm
param2xfm -translation -8 14 0 tmp/move_-8_15.xfm
param2xfm -translation -24 14 0 tmp/move_-24_15.xfm

mincresample tmp/rect_indent.mnc tmp/rect_indent_-8_15.mnc -transform tmp/move_-8_15.xfm -nearest -use_input_sampling
mincresample tmp/rect_indent.mnc tmp/rect_indent_24_15.mnc -transform tmp/move_24_15.xfm -nearest -use_input_sampling
mincresample tmp/rect_indent.mnc tmp/rect_indent_-24_15.mnc -transform tmp/move_-24_15.xfm -nearest -use_input_sampling
mincmath -max tmp/rect_indent_* tmp/indents_plus.mnc $calc
param2xfm -translation 0 -28 0 tmp/move_down.xfm -clobber 

param2xfm -translation 30 8 0 tmp/move_31_8.xfm -clobber
param2xfm -translation -30 8 0  tmp/move_-31_8.xfm -clobber
param2xfm -translation -30 -8 0 tmp/move_-31_-8.xfm -clobber
param2xfm -translation 30 -8 0  tmp/move_31_-8.xfm -clobber
param2xfm -translation 16 0 0  tmp/move_16.xfm
param2xfm -translation -16 0 0 tmp/move_-16.xfm

mincresample tmp/indents_plus.mnc -transform tmp/move_down.xfm -nearest -use_input_sampling tmp/indents_minus.mnc -clobb
make_phantom -rectangle -xwidth 0.5 -ywidth 0.5 -zwidth 17 $dimensions tmp/rect_indent.mnc -center 0 0 0 -clobber
mincresample tmp/rect_indent.mnc tmp/rect_indent_-8_15.mnc  -transform tmp/move_-8_15.xfm -nearest -use_input_sampling -clob
mincresample tmp/rect_indent.mnc tmp/rect_indent_8_15.mnc   -transform tmp/move_8_15.xfm -nearest -use_input_sampling -clob
mincresample tmp/rect_indent.mnc tmp/rect_indent_24_15.mnc  -transform tmp/move_24_15.xfm -nearest -use_input_sampling -clob
mincresample tmp/rect_indent.mnc tmp/rect_indent_-24_15.mnc -transform tmp/move_-24_15.xfm -nearest -use_input_sampling -clob
mincmath -max tmp/rect_indent_* tmp/indents_plus.mnc -clobber $calc


mincresample tmp/indents_plus.mnc -transform tmp/move_down.xfm -nearest -use_input_sampling tmp/indents_minus.mnc -clob

mincresample tmp/rect_indent.mnc tmp/rect_indent_31_8.mnc   -transform tmp/move_31_8.xfm -nearest -use_input_sampling -clob
mincresample tmp/rect_indent.mnc tmp/rect_indent_31_-8.mnc  -transform tmp/move_31_-8.xfm -nearest -use_input_sampling -clob
mincresample tmp/rect_indent.mnc tmp/rect_indent_-31_-8.mnc -transform tmp/move_-31_-8.xfm -nearest -use_input_sampling -clob
mincresample tmp/rect_indent.mnc tmp/rect_indent_-31_8.mnc  -transform tmp/move_-31_8.xfm -nearest -use_input_sampling -clob

mincmath -max tmp/rect_indent_* tmp/indents_plus.mnc tmp/indents_minus.mnc tmp/all_indents.mnc -clobber $calc

mincmath -max tmp/all_indents.mnc tmp/rect_shell.mnc tmp/rect_carcas.mnc $calc

make_phantom -ellipse -xwidth 13 -ywidth 13 -zwidth 2000 $dimensions tmp/ellipse.mnc -center 0 0 0 -clobber
make_phantom -ellipse -xwidth 11 -ywidth 11 -zwidth 2000 $dimensions tmp/ellipse_inner.mnc -center 0 0 0 -clobber
mincmath -sub tmp/ellipse.mnc tmp/ellipse_inner.mnc tmp/cylinder.mnc $calc

make_phantom -rectangle -xwidth 20 -ywidth 20 -zwidth 16 $dimensions tmp/ellipse_cut.mnc   -center 0 0 0 -clobber
minccalc -express 'A[1]>0.5?A[0]:0' tmp/cylinder.mnc tmp/ellipse_cut.mnc tmp/center_cylinder.mnc -clob  $calc

mincresample tmp/center_cylinder.mnc tmp/cyl_16.mnc  -transform  tmp/move_16.xfm -nearest -use_input_sampling -clob
mincresample tmp/center_cylinder.mnc tmp/cyl_-16.mnc -transform tmp/move_-16.xfm -nearest -use_input_sampling -clob
mincmath -max tmp/center_cylinder.mnc tmp/cyl_-16.mnc tmp/cyl_16.mnc  tmp/cylinders.mnc -clob $calc

# param2xfm -translation 0 0 -2 tmp/move_cyl.xfm
# mincresample -nearest -use_input_sampling tmp/cylinders.mnc tmp/cylinders2.mnc -transform tmp/move_cyl.xfm

make_phantom -rectangle -xwidth 0.5 -ywidth 8 -zwidth 12.50 $dimensions tmp/wall1.mnc -center 0  10.5 -1.50 -clobber
make_phantom -rectangle -xwidth 0.5 -ywidth 8 -zwidth 12.50 $dimensions tmp/wall2.mnc -center 0 -10.5 -1.50 -clobber

mincmath -max tmp/rect_carcas.mnc tmp/cylinders.mnc tmp/wall1.mnc tmp/wall2.mnc tmp/lower_brick.mnc -clob $calc
make_phantom -ellipse -xwidth 9 -ywidth 9 -zwidth 10000 $dimensions tmp/stud.mnc -center 0 0 0 
make_phantom -ellipse -xwidth 7 -ywidth 7 -zwidth 10000 $dimensions tmp/stud_inner.mnc -center 0 0 0 
mincmath -sub tmp/stud.mnc tmp/stud_inner.mnc tmp/stud_walls.mnc $calc

#                                    v stud height
make_phantom -rectangle -width 10 10 5 $dimensions tmp/stud_mask.mnc -center 0 0 -11.5 -clob

minccalc -express 'A[1]>0.5?A[0]:0' tmp/stud_walls.mnc tmp/stud_mask.mnc tmp/stud_c.mnc  $calc

param2xfm -translation -8 8 0  tmp/move_-8_8.xfm
param2xfm -translation -8 -8 0 tmp/move_-8_-8.xfm
param2xfm -translation 8 -8 0  tmp/move_8_-8.xfm
mincresample -nearest -use_input_sampling tmp/stud_c.mnc -transform tmp/move_8_8.xfm   tmp/stud_8_8.mnc
mincresample -nearest -use_input_sampling tmp/stud_c.mnc -transform tmp/move_-8_8.xfm  tmp/stud_-8_8.mnc
mincresample -nearest -use_input_sampling tmp/stud_c.mnc -transform tmp/move_-8_-8.xfm tmp/stud_-8_-8.mnc
mincresample -nearest -use_input_sampling tmp/stud_c.mnc -transform tmp/move_8_-8.xfm  tmp/stud_8_-8.mnc
mincmath -max tmp/stud_-8_-8.mnc tmp/stud_-8_8.mnc tmp/stud_8_-8.mnc tmp/stud_8_8.mnc tmp/studs.mnc $calc
mincresample tmp/studs.mnc -nearest tmp/studs_right.mnc -transform tmp/move_16.xfm -use_input_sampling
mincresample tmp/studs.mnc -nearest tmp/studs_left.mnc -transform  tmp/move_-16.xfm -use_input_sampling
mincmath -max tmp/studs_right.mnc tmp/studs_left.mnc tmp/studs_all.mnc $calc
make_phantom -rectangle -xwidth 3 -ywidth 0.5 -zwidth 12 tmp/wall_small1.mnc -center 8 0 -2 -clobber $dimensions
make_phantom -rectangle -xwidth 3 -ywidth 0.5 -zwidth 12 tmp/wall_small2.mnc -center -8 0 -2 -clobber $dimensions
mincmath -max tmp/studs_all.mnc tmp/lower_brick.mnc tmp/wall_small1.mnc tmp/wall_small2.mnc tmp/lego_brick.mnc $calc


export MINC_FORCE_V2=1
export MINC_COMPRESS=4

minccalc -expression 'clamp(A[0],0,1)' -byte tmp/lego_brick.mnc lego_brick_ideal.mnc      -clob $calc
minccalc -expression 'clamp(A[0],0,1)' -byte tmp/rect_outer.mnc lego_brick_ideal_mask.mnc -clob $calc
#param2xfm -translation 0 0 18 move_18.xfm 
#mincresample -transform move_18.xfm -use_input_sampling lego_brick_ideal.mnc lego_brick_ideal_moved.mnc -nearest

