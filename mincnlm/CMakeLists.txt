if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS} -DUSE_OPENMP")
    set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

ADD_EXECUTABLE( mincnlm
                mincnlm.cpp 
                nl_means.cpp
                nl_means_block.cpp
                nl_means_utils.cpp
                nl_means.h
                nl_means_block.h
                nl_means_utils.h )

TARGET_LINK_LIBRARIES( mincnlm minc_image_proc ${EZMINC_LIBRARIES} pthread  )

ADD_EXECUTABLE( minc_anlm
                minc_anlm.cpp 
              )
TARGET_LINK_LIBRARIES( minc_anlm minc_image_proc ${EZMINC_LIBRARIES} pthread  )

INSTALL(TARGETS 
        mincnlm minc_anlm
        RUNTIME DESTINATION bin )

