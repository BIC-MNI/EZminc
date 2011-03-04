ADD_LIBRARY(spherical_harmonics
  mincMeanSquaresImageToImageMetric.h
  data_proc.h
  sphericalHarmonicsTransform.cpp 
  sphericalHarmonicsTransform.h)

LINK_LIBRARIES(
    gsl
    gslcblas
    )

ADD_EXECUTABLE( fit_harmonics_grids
                fit_harmonics_grids.cpp )
                
TARGET_LINK_LIBRARIES(fit_harmonics_grids 
  minc4itk
  minc_io 
  spherical_harmonics gsl gslcblas)
  
ADD_EXECUTABLE( c_fit_harmonics_grids
                c_fit_harmonics_grids.cpp )
                
TARGET_LINK_LIBRARIES(c_fit_harmonics_grids 
  minc4itk
  minc_io 
  spherical_harmonics 
  gsl 
  gslcblas)
  
ADD_EXECUTABLE( param2grid
                param2grid.cpp )
                
                
TARGET_LINK_LIBRARIES( param2grid 
  minc4itk
  minc_io 
  spherical_harmonics gsl gslcblas)   

ADD_EXECUTABLE( c_param2grid
                c_param2grid.cpp )
                
TARGET_LINK_LIBRARIES( c_param2grid 
  minc4itk
  minc_io 
  spherical_harmonics gsl gslcblas)   
  
                

INSTALL(TARGETS 
    fit_harmonics_grids
    c_fit_harmonics_grids
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
