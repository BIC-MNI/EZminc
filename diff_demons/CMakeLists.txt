# Add include directories needed to use ITK.
INCLUDE_DIRECTORIES(BEFORE ${ITK_INCLUDE_DIRS})

# Add link directories needed to use ITK.
LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})

INCLUDE_DIRECTORIES( itk )

ADD_EXECUTABLE(DemonsRegistration DemonsRegistration.cxx mincUtils.cxx)
TARGET_LINK_LIBRARIES( DemonsRegistration ${VOLUME_IO_LIBRARIES} ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

ADD_EXECUTABLE(LogDomainDemonsRegistration LogDomainDemonsRegistration.cxx mincUtils.cxx)
TARGET_LINK_LIBRARIES( LogDomainDemonsRegistration ${VOLUME_IO_LIBRARIES} ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

ADD_EXECUTABLE(grid_2_log grid_2_log.cxx mincUtils.cxx)
TARGET_LINK_LIBRARIES( grid_2_log  ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

ADD_EXECUTABLE(log_resample log_resample.cpp)
TARGET_LINK_LIBRARIES( log_resample  ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

ADD_EXECUTABLE(grid_statistics grid_statistics.cpp)
TARGET_LINK_LIBRARIES( grid_statistics  ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

INSTALL(TARGETS DemonsRegistration LogDomainDemonsRegistration grid_2_log  log_resample grid_statistics RUNTIME DESTINATION bin )

