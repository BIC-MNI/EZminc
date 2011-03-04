/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: gaussian blurring using FFTW3
@COPYRIGHT  :
              Copyright 2009 Vladimir Fonov, 
              McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#ifndef __FFTW_BLUR_H__
#define __FFTW_BLUR_H__

#include <vector>
#include <limits>
#include <complex>
#include <minc_io_simple_volume.h>
namespace minc
{
  typedef std::complex<float> complex_f;
  //! calculate a windowed [derivative of ] gaussian
  void calculate_gaussian(std::vector<complex_f>& res,double sigma,bool dx);
  
  //! coordinate wise blurring
  void blur_volume(simple_volume<float> &in,simple_volume<float> &out, bool dx,bool dy, bool dz, double fx,double fy,double fz);
  
  void calc_gradient(simple_volume<float> &in,minc_grid_volume &out, double fx,double fy,double fz);
  void calc_gradient_mag(simple_volume<float> &in,simple_volume<float> &out, double fx,double fy,double fz);

};

#endif //__FFTW_BLUR_H__
