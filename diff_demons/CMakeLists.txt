# Add include directories needed to use ITK.
INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS} ${MINC4ITK_INCLUDE_DIRS})  
# Add link directories needed to use ITK.
LINK_DIRECTORIES(${ITK_LIBRARY_DIRS} ${MINC4ITK_LIBRARY_DIRS} ${MINC4ITK_INCLUDE_DIRS})

INCLUDE_DIRECTORIES( itk )

LINK_LIBRARIES( ITKIOMINC  ${ITK_LIBRARIES})

ADD_EXECUTABLE(DemonsRegistration DemonsRegistration.cxx mincUtils.cxx)

ADD_EXECUTABLE(LogDomainDemonsRegistration LogDomainDemonsRegistration.cxx mincUtils.cxx)

ADD_EXECUTABLE(grid_2_log grid_2_log.cxx mincUtils.cxx)

ADD_EXECUTABLE(log_resample log_resample.cpp)

ADD_EXECUTABLE(grid_statistics grid_statistics.cpp)

INSTALL(TARGETS DemonsRegistration 
                LogDomainDemonsRegistration 
                grid_2_log  
                log_resample 
                grid_statistics 
      RUNTIME DESTINATION bin )

