ADD_EXECUTABLE( mincnlm
                mincnlm.cpp 
                nl_means.cpp
                nl_means_block.cpp
                nl_means_utils.cpp
                nl_means.h
                nl_means_block.h
                nl_means_utils.h )

TARGET_LINK_LIBRARIES( mincnlm minc_image_proc minc_io pthread  )

INSTALL(TARGETS 
        mincnlm
        RUNTIME DESTINATION bin )

