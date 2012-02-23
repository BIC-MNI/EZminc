/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
              This routine implements simple algorithm of the sharpness estimate
              based on the median value of the gradient magnitude.

@COPYRIGHT  :
              Copyright 2012 Vladimir Fonov, 
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
#ifndef __SHARPNESS_ESTIMATE_H__
#define __SHARPNESS_ESTIMATE_H__

namespace minc
{
  double sharpness_estimate(const minc::simple_volume<float>& input,const minc::simple_volume<unsigned char>& mask);

};

#endif //__NOISE_ESTIMATE_H__