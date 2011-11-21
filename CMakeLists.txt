PROJECT (ezminc)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
 cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

INCLUDE_DIRECTORIES(
  minc_io image_proc
  )

SET(MINC_ROOT "/usr/local/bic" CACHE PATH "minc library base")
SET(NETCDF_ROOT "/usr/local" CACHE PATH "Netcdf library base")

OPTION(USE_MINC2      "Use minc 2 library" ON)
OPTION(BUILD_TESTS    "Build tests" OFF)
OPTION(BUILD_EXAMPLES "Build Examples" OFF)
OPTION(BUILD_MINCNLM  "Build Non-local means denoising programm,requires GSL (protected by a license IDDN.FR.001.070033.000.S.P.2007.000.21000)" OFF)

OPTION(BUILD_DISTORTION_CORRECTION "Build distortion correction tools, requires GSL" OFF)
OPTION(HAVE_GSL       "GSL is installed" OFF)
OPTION(HAVE_FFTW3     "FFTW3 is installed" OFF)

LINK_DIRECTORIES(
"${MINC_ROOT}/lib"
"${NETCDF_ROOT}/lib"
)
INCLUDE_DIRECTORIES( 
"${MINC_ROOT}/include"
"${NETCDF_ROOT}/include"
)


IF(HAVE_GSL)
  SET(GSL_ROOT "/usr" CACHE PATH "GSL library base")
  LINK_DIRECTORIES(
  "${GSL_ROOT}/lib"
  )

  INCLUDE_DIRECTORIES( 
    "${GSL_ROOT}/include"
    tools
  )

  IF(BUILD_MINCNLM)
    SUBDIRS(mincnlm)
  ENDIF(BUILD_MINCNLM)

ENDIF(HAVE_GSL)

IF(HAVE_FFTW3)

  SET(FFTW3_ROOT "/usr" CACHE PATH "FFTW3 library base")
  
  LINK_DIRECTORIES(
  "${FFTW3_ROOT}/lib"
  )
  INCLUDE_DIRECTORIES( 
  "${FFTW3_ROOT}/include"
  )
ENDIF(HAVE_FFTW3)


SUBDIRS( minc_io image_proc )

FIND_PACKAGE(ITK)

IF(ITK_FOUND)
  INCLUDE_DIRECTORIES(
    minc4itk
    ${ITK_INCLUDE_DIRS}
    )
  LINK_DIRECTORIES( ${ITK_LIBRARY_DIRS} )
  SUBDIRS(minc4itk)
  
  OPTION(BUILD_TOOLS    "Build additional tools (defacing) requires GSL" OFF)
ENDIF(ITK_FOUND)

IF(BUILD_EXAMPLES)
  SUBDIRS(examples)
ENDIF(BUILD_EXAMPLES)

IF(BUILD_TESTS)
  SUBDIRS(tests)
ENDIF(BUILD_TESTS)

IF(BUILD_TOOLS )
  SUBDIRS(tools)
ENDIF(BUILD_TOOLS)

IF(HAVE_GSL AND BUILD_DISTORTION_CORRECTION)
  SUBDIRS(distortion_correction)
ENDIF(HAVE_GSL AND BUILD_DISTORTION_CORRECTION)
