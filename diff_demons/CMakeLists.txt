INCLUDE_DIRECTORIES( itk )

ADD_EXECUTABLE(DemonsRegistration DemonsRegistration.cxx mincUtils.cxx)
TARGET_LINK_LIBRARIES( DemonsRegistration minc_tools ${VOLUME_IO_LIBRARIES} ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

ADD_EXECUTABLE(LogDomainDemonsRegistration LogDomainDemonsRegistration.cxx mincUtils.cxx)
TARGET_LINK_LIBRARIES( LogDomainDemonsRegistration minc_tools ${VOLUME_IO_LIBRARIES} ${ITK_LIBRARIES} ${MINC4ITK_LIBRARIES})

ADD_EXECUTABLE(grid_2_log grid_2_log.cxx mincUtils.cxx)
TARGET_LINK_LIBRARIES( grid_2_log  minc_tools)

ADD_EXECUTABLE(log_resample log_resample.cpp)
TARGET_LINK_LIBRARIES( log_resample  minc_tools)

ADD_EXECUTABLE(grid_statistics grid_statistics.cpp)
TARGET_LINK_LIBRARIES( grid_statistics  minc_tools)

INSTALL(TARGETS DemonsRegistration LogDomainDemonsRegistration grid_2_log  log_resample grid_statistics RUNTIME DESTINATION bin )


