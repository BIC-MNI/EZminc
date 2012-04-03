IF(EXISTS ${ITK_USE_FILE})
  INCLUDE(${ITK_USE_FILE})
ENDIF(EXISTS ${ITK_USE_FILE})

LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})
INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS} ${MINC4ITK_INCLUDE_DIRS})  

ADD_LIBRARY(spherical_harmonics
  mincMeanSquaresImageToImageMetric.h
  data_proc.h
  sphericalHarmonicsTransform.cpp 
  sphericalHarmonicsTransform.h)

TARGET_LINK_LIBRARIES( spherical_harmonics
    ITKIOMINC 
    ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} 
    )

ADD_EXECUTABLE( fit_harmonics_grids
                fit_harmonics_grids.cpp )
                
TARGET_LINK_LIBRARIES(fit_harmonics_grids 
  ITKIOMINC ${ITK_LIBRARIES} 
  spherical_harmonics )
  
ADD_EXECUTABLE( fit_harmonics_grids_regularize
                fit_harmonics_grids_regularize.cpp )
                
TARGET_LINK_LIBRARIES(fit_harmonics_grids_regularize 
  ITKIOMINC ${ITK_LIBRARIES}
  spherical_harmonics )
  
  
ADD_EXECUTABLE( c_fit_harmonics_grids
                c_fit_harmonics_grids.cpp )
                
TARGET_LINK_LIBRARIES(c_fit_harmonics_grids 
  spherical_harmonics 
  ${ITK_LIBRARIES}
  )
  
ADD_EXECUTABLE( param2grid
                param2grid.cpp )
                
                
TARGET_LINK_LIBRARIES( param2grid 
  spherical_harmonics 
  ${ITK_LIBRARIES}
  )

ADD_EXECUTABLE( sph_par2grid
                sph_par2grid.cpp )
                
                
TARGET_LINK_LIBRARIES( sph_par2grid 
  spherical_harmonics 
  ${ITK_LIBRARIES}
  )



ADD_EXECUTABLE( c_param2grid
                c_param2grid.cpp )
                
TARGET_LINK_LIBRARIES(c_param2grid 
  spherical_harmonics 
  ${ITK_LIBRARIES}
  )
                  

INSTALL(TARGETS 
    fit_harmonics_grids
    c_fit_harmonics_grids
    fit_harmonics_grids_regularize
    param2grid
    c_param2grid
    RUNTIME DESTINATION bin )

INSTALL(PROGRAMS 
    calc_distortions.pl
    par2xfm.pl
    phantom_distortion_measure.pl
    lego_core_extract.pl
    adni_preprocess.pl    
    phantomfit.pl
    DESTINATION bin )
