#ifndef _MINC_WRAPPERS_H_
#define _MINC_WRAPPERS_H_

#include <complex>
#include <vector>
#include <algorithm>
#include <itkArray.h>
#include <iostream>

#include <minc_io_exceptions.h>
#include <minc_io_fixed_vector.h>

#if (ITK_VERSION_MAJOR < 4)
#include <itkMincHelpers.h>
#else
#include "itk4MincHelpers.h"
#endif

#endif //_MINC_WRAPPERS_H_
