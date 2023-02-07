LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})

INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})


message("ITK_LIBRARIES=${ITK_LIBRARIES}")

IF( HAVE_MINC4ITK)
  INCLUDE_DIRECTORIES(${MINC4ITK_INCLUDE_DIRS})
ENDIF()

ADD_LIBRARY(spherical_harmonics
  mincMeanSquaresImageToImageMetric.h
  data_proc.h
  sphericalHarmonicsTransform.cpp
  sphericalHarmonicsTransform.h)

  IF( ITK_VERSION_MAJOR VERSION_GREATER_EQUAL 4 )
  #SET(ITK_LIBRARIES itk4MincHelpers
  SET(ITK_LIBRARIES itk4MincHelpers ${ITK_LIBRARIES})
  ENDIF( ITK_VERSION_MAJOR VERSION_GREATER_EQUAL 4 )

  TARGET_LINK_LIBRARIES(
    spherical_harmonics
    minc_image_proc

    ${ITK_LIBRARIES}
    ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} 
    )

ADD_EXECUTABLE( fit_harmonics_grids
                fit_harmonics_grids.cpp )

TARGET_LINK_LIBRARIES(fit_harmonics_grids 
  spherical_harmonics ${ITK_LIBRARIES} ${ZLIB_LIBRARY} )

ADD_EXECUTABLE( fit_harmonics_grids_diff
                fit_harmonics_grids_diff.cpp )

TARGET_LINK_LIBRARIES(fit_harmonics_grids_diff
  spherical_harmonics
  ${ITK_LIBRARIES} ${ZLIB_LIBRARY} )


ADD_EXECUTABLE( fit_harmonics_grids_regularize
                fit_harmonics_grids_regularize.cpp )

TARGET_LINK_LIBRARIES(fit_harmonics_grids_regularize 
  spherical_harmonics
  ${ITK_LIBRARIES} ${ZLIB_LIBRARY}
  )

ADD_EXECUTABLE( c_fit_harmonics_grids
                c_fit_harmonics_grids.cpp )

TARGET_LINK_LIBRARIES(c_fit_harmonics_grids 
  spherical_harmonics
  ${ITK_LIBRARIES} ${ZLIB_LIBRARY}
  )

ADD_EXECUTABLE( param2grid
                param2grid.cpp )

TARGET_LINK_LIBRARIES( param2grid 
  spherical_harmonics
  ${ITK_LIBRARIES} ${ZLIB_LIBRARY}
  )

ADD_EXECUTABLE( sph_par2grid
                sph_par2grid.cpp )
                
                
TARGET_LINK_LIBRARIES( sph_par2grid 
  spherical_harmonics 
  ${ITK_LIBRARIES} ${ZLIB_LIBRARY}
  )



ADD_EXECUTABLE( c_param2grid
                c_param2grid.cpp )
                
TARGET_LINK_LIBRARIES(c_param2grid 
  spherical_harmonics 
  ${ITK_LIBRARIES} ${ZLIB_LIBRARY}
  )

ADD_EXECUTABLE( itk_label_stats
                itk_label_stats.cpp )
                
TARGET_LINK_LIBRARIES( itk_label_stats 
  ${ITK_LIBRARIES} ${ZLIB_LIBRARY}
  )


INSTALL(TARGETS
    fit_harmonics_grids
    fit_harmonics_grids_diff
    c_fit_harmonics_grids
    fit_harmonics_grids_regularize
    param2grid
    c_param2grid
    itk_label_stats
#    itk_assemble_bricks
    RUNTIME DESTINATION bin )

INSTALL(PROGRAMS 
    calc_distortions.pl
    par2xfm.pl
    phantom_distortion_measure.pl
    phantom_distortion_measure_v2.pl
    lego_core_extract.pl
    adni_preprocess.pl    
    phantomfit.pl
    phantomfit_ANTS.pl
    phantomfit_DD.pl
    phantomfit_elastix.pl
    DESTINATION bin )
