ADD_LIBRARY( mrflib 
             generic_image.cpp  generic_image.h
             minc_generic_image.cpp minc_generic_image.h
             atlasspec.cpp  atlasspec.h  
             nhmrf.cpp      nhmrf.h 
             parzen.cpp     parzen.h )


ADD_EXECUTABLE( mrfseg
                mrfseg.cpp )
                
TARGET_LINK_LIBRARIES( mrfseg  
                mrflib ${EZMINC_LIBRARIES})

ADD_EXECUTABLE( gamixture
                gamixture.cpp gamixture.h 
                gamixturemain.cpp)
                
TARGET_LINK_LIBRARIES(gamixture 
                      mrflib ${EZMINC_LIBRARIES})
                
ADD_EXECUTABLE(relabel_analyze
               relabel_analyze.cpp)
                 
TARGET_LINK_LIBRARIES(relabel_analyze
                      mrflib ${EZMINC_LIBRARIES})


INSTALL(TARGETS 
#        relabel_analyze
        gamixture 
        mrfseg
      RUNTIME DESTINATION bin )

INSTALL(PROGRAMS 
  mrfseg.pl
  DESTINATION bin )
