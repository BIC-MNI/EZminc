IF(ITK_FOUND)

  ADD_EXECUTABLE( minc_io_test_anat 
                  minc_io_test_anat.cpp)
  TARGET_LINK_LIBRARIES( minc_io_test_anat 
                minc4itk)
                
  ADD_EXECUTABLE( any_itk_any 
                  any_itk_any.cpp)
  TARGET_LINK_LIBRARIES( any_itk_any 
                minc4itk)
                
  ADD_EXECUTABLE( any_itk_any_4d
                  any_itk_any_4d.cpp)
  TARGET_LINK_LIBRARIES( any_itk_any_4d 
                minc4itk)
                
  ADD_EXECUTABLE( any_itk_any_u32
                  any_itk_any_u32.cpp)
  TARGET_LINK_LIBRARIES( any_itk_any_u32
                minc4itk)
                
  ADD_EXECUTABLE( minc_io_test_grid 
                  minc_io_test_grid.cpp)
  TARGET_LINK_LIBRARIES( minc_io_test_grid
                minc4itk)
                
  ADD_EXECUTABLE( minc_io_test_dti 
                  minc_io_test_dti.cpp)
  TARGET_LINK_LIBRARIES( minc_io_test_dti  
                minc4itk)
  
  ADD_EXECUTABLE( minc_io_test_coordinates
                  minc_io_test_coordinates.cpp)
  TARGET_LINK_LIBRARIES( minc_io_test_coordinates  
                minc4itk)
                
  ADD_EXECUTABLE( minc_io_test_coordinates2
                  minc_io_test_coordinates2.cpp)
  TARGET_LINK_LIBRARIES( minc_io_test_coordinates2
                minc4itk)
                
  ADD_EXECUTABLE( tags_rw_test
                  tags_rw_test.cpp)
  TARGET_LINK_LIBRARIES( tags_rw_test
                minc4itk)
    
ENDIF(ITK_FOUND)

ADD_EXECUTABLE( minc_rw_test
                minc_rw_test.cpp )
TARGET_LINK_LIBRARIES( minc_rw_test minc_io)

ADD_EXECUTABLE( minc_rw_test2
                minc_rw_test2.cpp )
TARGET_LINK_LIBRARIES( minc_rw_test2 minc_io)


ADD_EXECUTABLE( minc_rw_test_simple
                minc_rw_test_simple.cpp )
TARGET_LINK_LIBRARIES( minc_rw_test_simple minc_io)

ADD_EXECUTABLE( minc_rw_test_4d
               minc_rw_test_4d.cpp )
TARGET_LINK_LIBRARIES( minc_rw_test_4d minc_io)


ADD_EXECUTABLE( transformpoint
                transformpoint.cpp )
TARGET_LINK_LIBRARIES( transformpoint minc_io)


ADD_EXECUTABLE( minc_com
                minc_com.cpp )
TARGET_LINK_LIBRARIES( minc_com minc_io)

ADD_EXECUTABLE( minc_com_4d
                minc_com_4d.cpp )
TARGET_LINK_LIBRARIES( minc_com_4d minc_io)
