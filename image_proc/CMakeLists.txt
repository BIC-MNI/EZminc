
SET( IMAGE_PROC_SRC minc_histograms.cpp dwt.cpp  dwt_utils.cpp )
SET( IMAGE_PROC_HDR minc_histograms.h dwt.h  dwt_utils.h )

IF( HAVE_GSL )
	SET( IMAGE_PROC_HDR ${IMAGE_PROC_HDR} gsl_gauss.h  gsl_glue.h gsl_nlfit.h noise_estimate.h)
	SET( IMAGE_PROC_SRC ${IMAGE_PROC_SRC} gsl_nlfit.cpp noise_estimate.cpp)
ENDIF( HAVE_GSL )


IF( HAVE_FFTW3 )
	SET( IMAGE_PROC_HDR ${IMAGE_PROC_HDR} fftw_blur.h)
	SET( IMAGE_PROC_SRC ${IMAGE_PROC_SRC} fftw_blur.cpp)
ENDIF( HAVE_FFTW3 )


ADD_LIBRARY( minc_image_proc ${IMAGE_PROC_HDR} ${IMAGE_PROC_SRC})

IF(HAVE_GSL)
	TARGET_LINK_LIBRARIES(minc_image_proc gsl gslcblas)
ENDIF(HAVE_GSL)

IF(HAVE_FFTW3)
	TARGET_LINK_LIBRARIES(minc_image_proc fftw3f)
ENDIF(HAVE_FFTW3)
