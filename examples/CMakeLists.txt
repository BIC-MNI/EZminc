IF(ITK_FOUND)

ADD_EXECUTABLE( volume_2_csv volume_2_csv.cpp)
TARGET_LINK_LIBRARIES( volume_2_csv  minc4itk)

ADD_EXECUTABLE( itk_resample 
                itk_resample.cpp )
TARGET_LINK_LIBRARIES( itk_resample minc4itk)

ADD_EXECUTABLE( volume_moments 
                volume_moments.cpp )
TARGET_LINK_LIBRARIES( volume_moments minc4itk)

ADD_EXECUTABLE( itk_dti 
                itk_dti.cpp )
TARGET_LINK_LIBRARIES( itk_dti minc4itk)

INSTALL(TARGETS 
        itk_resample volume_moments
        RUNTIME DESTINATION bin )
        
ENDIF(ITK_FOUND)

ADD_EXECUTABLE( volume_msq_dist volume_msq_dist.cpp)
TARGET_LINK_LIBRARIES( volume_msq_dist minc_io)

ADD_EXECUTABLE( volume_avg volume_avg.cpp)
TARGET_LINK_LIBRARIES( volume_avg minc_io)

ADD_EXECUTABLE( volume_similarity volume_similarity.cpp)
TARGET_LINK_LIBRARIES( volume_similarity minc_io)

ADD_EXECUTABLE( trilinear_resample trilinear_resample.cpp)
TARGET_LINK_LIBRARIES( trilinear_resample minc_io)


INSTALL(TARGETS 
        volume_msq_dist
        volume_avg
        volume_similarity
        RUNTIME DESTINATION bin )
        
