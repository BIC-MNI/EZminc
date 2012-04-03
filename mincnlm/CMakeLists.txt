ADD_EXECUTABLE( mincnlm
                mincnlm.cpp 
                nl_means.cpp
                nl_means_block.cpp
                nl_means_utils.cpp
                nl_means.h
                nl_means_block.h
                nl_means_utils.h )

TARGET_LINK_LIBRARIES( mincnlm minc_image_proc ${EZMINC_LIBRARIES} pthread  )

# IF(GSL_FOUND)
#   TARGET_LINK_LIBRARIES(mincnlm ${GSL_LIBRARIES})
# ENDIF(GSL_FOUND)
# 
# IF(FFTW3F_FOUND)
#   TARGET_LINK_LIBRARIES(mincnlm ${FFTW3F_LIBRARIES})
# ENDIF(FFTW3F_FOUND)



INSTALL(TARGETS 
        mincnlm
        RUNTIME DESTINATION bin )

