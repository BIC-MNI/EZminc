IF(GSL_FOUND)
	ADD_EXECUTABLE(volume_pol
			volume_pol.cpp )

	TARGET_LINK_LIBRARIES( volume_pol 
		minc_image_proc ${EZMINC_LIBRARIES}  gsl gslcblas)
		
	
  ADD_EXECUTABLE( random_volume
									random_volume.cpp )

	TARGET_LINK_LIBRARIES(random_volume ${EZMINC_LIBRARIES}  gsl gslcblas)

	INSTALL(TARGETS 
        random_volume
        volume_pol
        RUNTIME DESTINATION bin )
ENDIF(GSL_FOUND)


IF(FFTW3F_FOUND)
	ADD_EXECUTABLE( fast_blur
			fast_blur.cpp )

	TARGET_LINK_LIBRARIES( fast_blur
		minc_image_proc ${EZMINC_LIBRARIES})
		
	INSTALL(TARGETS 
        fast_blur
        RUNTIME DESTINATION bin )
ENDIF(FFTW3F_FOUND)


IF(FFTW3F_FOUND AND GSL_FOUND)
  ADD_EXECUTABLE( noise_estimate noise_estimate.cpp)

  TARGET_LINK_LIBRARIES(noise_estimate minc_image_proc ${EZMINC_LIBRARIES})

  ADD_EXECUTABLE( sharpness_estimate sharpness_estimate.cpp)

  TARGET_LINK_LIBRARIES(sharpness_estimate minc_image_proc ${EZMINC_LIBRARIES})
  
  
  INSTALL(TARGETS 
        noise_estimate sharpness_estimate
        RUNTIME DESTINATION bin )

  INSTALL(PROGRAMS 
        volume_denoise.pl
        DESTINATION bin )
ENDIF(FFTW3F_FOUND AND GSL_FOUND)


IF(ITK_FOUND AND HAVE_MINC4ITK)
  INCLUDE(${ITK_USE_FILE})

  LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})
  INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})  

  ADD_EXECUTABLE( resample_grid
                resample_grid.cpp )

  TARGET_LINK_LIBRARIES(resample_grid ${ITK_LIBRARIES} ${VOLUME_IO_LIBRARIES} ${MINC4ITK_LIBRARIES})
	
  ADD_EXECUTABLE( multiple_volume_similarity multiple_volume_similarity.cpp)
  TARGET_LINK_LIBRARIES( multiple_volume_similarity ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES} )
        
  ADD_EXECUTABLE( tag2csv  tag2csv.cpp)
  TARGET_LINK_LIBRARIES(tag2csv ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_distance itk_distance.cpp)
	
  INSTALL(TARGETS 
   resample_grid multiple_volume_similarity itk_distance
   RUNTIME DESTINATION bin )

ENDIF(ITK_FOUND AND HAVE_MINC4ITK)


ADD_EXECUTABLE( volume_dwt volume_dwt.cpp)

TARGET_LINK_LIBRARIES(volume_dwt minc_image_proc ${EZMINC_LIBRARIES})

ADD_EXECUTABLE( minc_downsample 
                minc_downsample.cpp )
TARGET_LINK_LIBRARIES(minc_downsample ${EZMINC_LIBRARIES})

ADD_EXECUTABLE( minc_nuyl
                minc_nuyl.cpp )
TARGET_LINK_LIBRARIES(minc_nuyl minc_image_proc ${EZMINC_LIBRARIES})


ADD_EXECUTABLE( minc_rank
                minc_rank.cpp )
TARGET_LINK_LIBRARIES(minc_rank minc_image_proc ${EZMINC_LIBRARIES})


INSTALL(TARGETS 
        volume_dwt minc_downsample minc_nuyl minc_rank
        RUNTIME DESTINATION bin )

INSTALL(PROGRAMS 
        deface_minipipe.pl
        make_random_grid.pl
        deface_volume.pl
        DESTINATION bin )
