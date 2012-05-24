IF(ITK_FOUND)
  INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS} ${MINC4ITK_INCLUDE_DIRS})  
  LINK_DIRECTORIES(${ITK_LIBRARY_DIRS} ${MINC4ITK_LIBRARY_DIRS})
  
  ADD_EXECUTABLE( volume_2_csv volume_2_csv.cpp)
  TARGET_LINK_LIBRARIES( volume_2_csv  ITKIOMINC ${ITK_LIBRARIES})

  ADD_EXECUTABLE( volume_moments 
                  volume_moments.cpp )
  TARGET_LINK_LIBRARIES( volume_moments ITKIOMINC ${ITK_LIBRARIES})

  ADD_EXECUTABLE( itk_dti 
                  itk_dti.cpp )
  TARGET_LINK_LIBRARIES( itk_dti ITKIOMINC ${ITK_LIBRARIES})
  
  ADD_EXECUTABLE( itk_test
                  itk_test.cpp )
  TARGET_LINK_LIBRARIES( itk_test ITKIOMINC ${ITK_LIBRARIES})

  INSTALL(TARGETS 
          volume_moments
          RUNTIME DESTINATION bin )
  
  IF(BICPL_FOUND)
    INCLUDE_DIRECTORIES(${BICPL_INCLUDE_DIRS})  
    ADD_EXECUTABLE( object_intersect object_intersect.cpp )
    TARGET_LINK_LIBRARIES( object_intersect ${BICPL_LIBRARIES} ITKIOMINC ${ITK_LIBRARIES})
  ENDIF(BICPL_FOUND)
        
ENDIF(ITK_FOUND)
