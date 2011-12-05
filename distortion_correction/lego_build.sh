#! /bin/sh
make_phantom -rectangle -xwidth 62 -ywidth 30 -zwidth 18 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 rect_outer.mnc -center 0 0 0
make_phantom -rectangle -xwidth 60 -ywidth 28 -zwidth 16 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 rect_inner.mnc -center 0 0 1
mincmath -sub rect_outer.mnc rect_inner.mnc rect_shell.mnc

make_phantom -rectangle -xwidth 1 -ywidth 1 -zwidth 17 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 rect_indent.mnc -center 0 0 0
param2xfm -translation 8 8 0 move_8_8.xfm
mincresample rect_indent.mnc rect_indent_8_9.mnc -transform move_8_8.xfm -nearest -use_input_sampling
rm rect_indent_8_9.mnc move_8_8.xfm
param2xfm -translation 14 8 0 move_15_8.xfm
mincresample rect_indent.mnc rect_indent_15_8.mnc -transform move_15_8.xfm -nearest -use_input_sampling
rm rect_indent_15_8.mnc move_15_8.xfm 
param2xfm -translation 8 14 0 move_8_15.xfm
mincresample rect_indent.mnc rect_indent_8_15.mnc -transform move_8_15.xfm -nearest -use_input_sampling
param2xfm -translation 24 14 0 move_24_15.xfm
param2xfm -translation -8 14 0 move_-8_15.xfm
param2xfm -translation -24 14 0 move_-24_15.xfm
mincresample rect_indent.mnc rect_indent_-8_15.mnc -transform move_-8_15.xfm -nearest -use_input_sampling
mincresample rect_indent.mnc rect_indent_24_15.mnc -transform move_24_15.xfm -nearest -use_input_sampling
mincresample rect_indent.mnc rect_indent_-24_15.mnc -transform move_-24_15.xfm -nearest -use_input_sampling
mincmath -max rect_indent_* indents_plus.mnc
param2xfm -translation 0 -28 0 move_down.xfm -clobber
mincresample indents_plus.mnc -transform move_down.xfm -nearest -use_input_sampling indents_minus.mnc -clobb
make_phantom -rectangle -xwidth 0.5 -ywidth 0.5 -zwidth 17 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 rect_indent.mnc -center 0 0 0 -clobber
mincresample rect_indent.mnc rect_indent_-8_15.mnc -transform move_-8_15.xfm -nearest -use_input_sampling -clob
mincresample rect_indent.mnc rect_indent_8_15.mnc -transform move_8_15.xfm -nearest -use_input_sampling -clob
mincresample rect_indent.mnc rect_indent_24_15.mnc -transform move_24_15.xfm -nearest -use_input_sampling -clob
mincresample rect_indent.mnc rect_indent_-24_15.mnc -transform move_-24_15.xfm -nearest -use_input_sampling -clob
mincmath -max rect_indent_* indents_plus.mnc -clobber
mincresample indents_plus.mnc -transform move_down.xfm -nearest -use_input_sampling indents_minus.mnc -clob
param2xfm -translation 0 -28 0 move_down.xfm -clobber
mincresample indents_plus.mnc -transform move_down.xfm -nearest -use_input_sampling indents_minus.mnc -clob
param2xfm -translation 30 8 0 move_31_8.xfm -clobber
mincresample rect_indent.mnc rect_indent_31_8.mnc -transform move_31_8.xfm -nearest -use_input_sampling -clob
param2xfm -translation -30 8 0 move_-31_8.xfm -clobber
param2xfm -translation -30 -8 0 move_-31_-8.xfm -clobber
param2xfm -translation 30 -8 0 move_31_-8.xfm -clobber
mincresample rect_indent.mnc rect_indent_31_8.mnc -transform move_31_8.xfm -nearest -use_input_sampling -clob
mincresample rect_indent.mnc rect_indent_31_-8.mnc -transform move_31_-8.xfm -nearest -use_input_sampling -clob
mincresample rect_indent.mnc rect_indent_-31_-8.mnc -transform move_-31_-8.xfm -nearest -use_input_sampling -clob
mincresample rect_indent.mnc rect_indent_-31_8.mnc -transform move_-31_8.xfm -nearest -use_input_sampling -clob
mincmath -max rect_indent_* indents_plus.mnc indents_minus.mnc all_indents.mnc -clobber
mincmath -max all_indents.mnc rect_shell.mnc rect_carcas.mnc

make_phantom -ellipse -xwidth 8 -ywidth 8 -zwidth 17 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 ellipse.mnc -center 0 0 0 -clobber
make_phantom -ellipse -xwidth 13 -ywidth 8 -zwidth 17 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 ellipse.mnc -center 0 0 0 -clobber
make_phantom -ellipse -xwidth 13 -ywidth 13 -zwidth 17 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 ellipse.mnc -center 0 0 0 -clobber
make_phantom -ellipse -xwidth 13 -ywidth 13 -zwidth 2000 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 ellipse.mnc -center 0 0 0 -clobber
make_phantom -ellipse -xwidth 11 -ywidth 11 -zwidth 2000 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 ellipse_inner.mnc -center 0 0 0 -clobber
mincmath -sub ellipse.mnc ellipse_inner.mnc cylinder.mnc
minccalc -express 'A[1]>0.5?A[0]:0' cylinder.mnc rect_outer.mnc center_cylinder.mnc
param2xfm -translation 16 0 0 move_16.xfm
param2xfm -translation -16 0 0 move_-16.xfm
mincresample cylinder.mnc cyl_16.mnc -transform move_16.xfm -nearest -use_input_sampling -clob
mincresample cylinder.mnc cyl_-16.mnc -transform move_-16.xfm -nearest -use_input_sampling -clob
mincmath -max rect_carcas.mnc cylinder.mnc cyl_16.mnc cyl_-16.mnc carcas2.mnc
minccalc -express 'A[1]>0.5?A[0]:0' cylinder.mnc rect_inner.mnc center_cylinder.mnc -clob
minccalc -expression 'A[0]/0.88' center_cylinder.mnc center_cylinder1.mnc
mv center_cylinder1.mnc center_cylinder.mnc
mincresample center_cylinder.mnc cyl_-16.mnc -transform move_-16.xfm -nearest -use_input_sampling -clob
mincresample center_cylinder.mnc cyl_16.mnc -transform move_16.xfm -nearest -use_input_sampling -clob
mincmath -max rect_carcas.mnc center_cylinder.mnc cyl_-16.mnc cyl_16.mnc carcas2.mnc -clob
make_phantom -rectangle -xwidth 20 -ywidth 20 -zwidth 15 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 ellipse_inner.mnc -center 0 0 1 -clobber
make_phantom -rectangle -xwidth 20 -ywidth 20 -zwidth 15 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 ellipse_cut.mnc -center 0 0 1 -clobber
minccalc -express 'A[1]>0.5?A[0]:0' cylinder.mnc ellipse_cut.mnc center_cylinder.mnc -clob
mincresample center_cylinder.mnc cyl_16.mnc -transform move_16.xfm -nearest -use_input_sampling -clob
mincresample center_cylinder.mnc cyl_-16.mnc -transform move_-16.xfm -nearest -use_input_sampling -clob
mincmath -max rect_carcas.mnc center_cylinder.mnc cyl_-16.mnc cyl_16.mnc carcas2.mnc -clob
mincmath -max center_cylinder.mnc cyl_-16.mnc cyl_16.mnc cylinders.mnc -clob
param2xfm -translation 0 0 -2 move_cyl.xfm
mincresample -nearest -use_input_sampling cylinders.mnc cylinders2.mnc -transform move_cyl.xfm
make_phantom -rectangle -xwidth 0.5 -ywidth 8 -zwidth 12 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 wall1.mnc -center 0  11 -1 -clobber
make_phantom -rectangle -xwidth 0.5 -ywidth 8 -zwidth 12 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 wall2.mnc -center 0 -11 -1 -clobber
mincmath -max rect_carcas.mnc cylinders2.mnc wall1.mnc wall2.mnc lower_brick.mnc -clob
make_phantom -ellipse -xwidth 9 -ywidth 9 -zwidth 10000 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 stud.mnc -center 0 0 0 
make_phantom -ellipse -xwidth 7 -ywidth 7 -zwidth 10000 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 stud_inner.mnc -center 0 0 0 
mincmath -sub stud.mnc stud_inner.mnc stud_walls.mnc

#                                    v stud height
make_phantom -rectangle -width 10 10 4 -start -50 -50 -50 -nelements 100 100 100  -step 1 1 1 stud_mask.mnc -center 0 0 -12.5 -clob

minccalc -express 'A[1]>0.5?A[0]:0' stud_walls.mnc stud_mask.mnc stud_c.mnc 

param2xfm -translation 8 8 0 move_8_8.xfm
param2xfm -translation -8 8 0 move_-8_8.xfm
param2xfm -translation -8 -8 0 move_-8_-8.xfm
param2xfm -translation 8 -8 0 move_8_-8.xfm
mincresample -nearest -use_input_sampling stud_c.mnc -transform move_8_8.xfm stud_8_8.mnc
mincresample -nearest -use_input_sampling stud_c.mnc -transform move_-8_8.xfm stud_-8_8.mnc
mincresample -nearest -use_input_sampling stud_c.mnc -transform move_-8_-8.xfm stud_-8_-8.mnc
mincresample -nearest -use_input_sampling stud_c.mnc -transform move_8_-8.xfm stud_8_-8.mnc
mincmath -max stud_-8_-8.mnc stud_-8_8.mnc stud_8_-8.mnc stud_8_8.mnc studs.mnc
mincresample studs.mnc -nearest studs_right.mnc -transform move_16.xfm -use_input_sampling
mincresample studs.mnc -nearest studs_left.mnc -transform move_-16.xfm -use_input_sampling
mincmath -max studs_right.mnc studs_left.mnc studs_all.mnc
make_phantom -rectangle -xwidth 3 -ywidth 0.5 -zwidth 12 -start -50 -50 -50 -nelements 100 100 100 -step 1 1 1 wall_small1.mnc -center 8 0 -2 -clobber
make_phantom -rectangle -xwidth 3 -ywidth 0.5 -zwidth 12 -start -50 -50 -50 -nelements 100 100 100 -step 1 1 1 wall_small2.mnc -center -8 0 -2 -clobber
mincmath -max studs_all.mnc lower_brick.mnc wall_small1.mnc wall_small2.mnc lego_brick.mnc
minccalc -expression 'clamp(A[0],0,1)' lego_brick.mnc lego_brick_ideal.mnc
param2xfm -translation 0 0 18 move_18.xfm 
#mincresample -transform move_18.xfm -use_input_sampling lego_brick_ideal.mnc lego_brick_ideal_moved.mnc -nearest

