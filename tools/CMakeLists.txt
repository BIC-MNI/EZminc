IF(HAVE_GSL AND ITK_FOUND)
	ADD_EXECUTABLE(volume_pol
			volume_pol.cpp )

	TARGET_LINK_LIBRARIES( volume_pol 
		minc4itk minc_image_proc)
		
	INSTALL(TARGETS 
        volume_pol
        RUNTIME DESTINATION bin )
        
ENDIF(HAVE_GSL AND ITK_FOUND)

IF(HAVE_GSL)
	ADD_EXECUTABLE( random_volume
									random_volume.cpp )

	TARGET_LINK_LIBRARIES(random_volume minc_io  gsl gslcblas)

	INSTALL(TARGETS 
        random_volume
        RUNTIME DESTINATION bin )
ENDIF(HAVE_GSL)


IF(HAVE_FFTW3)
	ADD_EXECUTABLE( fast_blur
			fast_blur.cpp )

	TARGET_LINK_LIBRARIES( fast_blur
		minc_image_proc minc_io)
		
	INSTALL(TARGETS 
        fast_blur
        RUNTIME DESTINATION bin )
ENDIF(HAVE_FFTW3)


IF(HAVE_FFTW3 AND HAVE_GSL)
  ADD_EXECUTABLE( noise_estimate noise_estimate.cpp)

  TARGET_LINK_LIBRARIES(noise_estimate minc_image_proc minc_io)

  INSTALL(TARGETS 
        noise_estimate
        RUNTIME DESTINATION bin )

  INSTALL(PROGRAMS 
        volume_denoise.pl
        DESTINATION bin )
ENDIF(HAVE_FFTW3 AND HAVE_GSL)


IF(ITK_FOUND)
  ADD_EXECUTABLE( resample_grid
                resample_grid.cpp )

	TARGET_LINK_LIBRARIES(resample_grid minc4itk)
	
  ADD_EXECUTABLE( itk_convert
                itk_convert.cpp )

	TARGET_LINK_LIBRARIES(itk_convert minc4itk)
	
	
	INSTALL(TARGETS 
        resample_grid itk_convert
        RUNTIME DESTINATION bin )
        
  ADD_EXECUTABLE( tag2csv  tag2csv.cpp)
  TARGET_LINK_LIBRARIES(tag2csv minc4itk)
ENDIF(ITK_FOUND)


ADD_EXECUTABLE( volume_dwt volume_dwt.cpp)

TARGET_LINK_LIBRARIES(volume_dwt minc_image_proc minc_io)


INSTALL(TARGETS 
        volume_dwt
        RUNTIME DESTINATION bin )

INSTALL(PROGRAMS 
        deface_minipipe.pl
        make_random_grid.pl
        deface_volume.pl
        DESTINATION bin )
