
SET( IMAGE_PROC_SRC minc_histograms.cpp dwt.cpp  dwt_utils.cpp )
SET( IMAGE_PROC_HDR minc_histograms.h dwt.h  dwt_utils.h )

IF( GSL_FOUND )
	SET( IMAGE_PROC_HDR ${IMAGE_PROC_HDR} gsl_gauss.h  gsl_glue.h gsl_nlfit.h noise_estimate.h)
	SET( IMAGE_PROC_SRC ${IMAGE_PROC_SRC} gsl_nlfit.cpp noise_estimate.cpp)
ENDIF( GSL_FOUND )


IF( FFTW3F_FOUND )
	SET( IMAGE_PROC_HDR ${IMAGE_PROC_HDR} fftw_blur.h)
	SET( IMAGE_PROC_SRC ${IMAGE_PROC_SRC} fftw_blur.cpp)
ENDIF( FFTW3F_FOUND )


ADD_LIBRARY( minc_image_proc ${IMAGE_PROC_HDR} ${IMAGE_PROC_SRC})

IF(GSL_FOUND)
	TARGET_LINK_LIBRARIES(minc_image_proc ${GSL_LIBRARIES})
ENDIF(GSL_FOUND)

IF(FFTW3F_FOUND)
	TARGET_LINK_LIBRARIES(minc_image_proc ${FFTW3F_LIBRARIES})
ENDIF(FFTW3F_FOUND)

TARGET_LINK_LIBRARIES(minc_image_proc ${EZMINC_LIBRARIES})
