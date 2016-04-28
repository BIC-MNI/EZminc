INSTALL(PROGRAMS 
        make_face.pl
        volume_denoise.pl
        deface_minipipe.pl
        make_random_grid.pl
        deface_volume.pl
        nlfit_s 
        nlfit_l
        nlfit_o2
        nlfit_f
        bestlinreg.pl 
        bestlinreg_s 
        bestlinreg_g
        bestlinreg_s2 
        icc_mask.pl
        uniformize_minc.pl
        xfm_normalize.pl
#        xfmavg
        lobe_segment.pl
        lobes_to_volumes.pl
        minc_qc2.pl
        minc_qc.pl
        minc_qc_rgb.pl
        minc_qc_t2t1.pl
        minc_pretty_pic.pl
        minc_pretty_pic_m.pl
        DESTINATION bin )

# make a copy of lobe_segment.pl
INSTALL(PROGRAMS
          lobe_segment.pl
        DESTINATION bin 
        RENAME 
          lobe_segment)
