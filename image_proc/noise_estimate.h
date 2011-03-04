/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
							This routine implements noise estimation algorithm published in 
							Pierrick Coupe, Jose V. Manjon, Elias Gedamu, Douglas L. Arnold,
							Montserrat Robles, D. Louis Collins: An Object-Based Method for Rician
							Noise Estimation in MR Images. MICCAI (1) 2009: 601-608.

@COPYRIGHT  :
              Copyright 2009 Pierrick Coupe,Vladimir Fonov, 
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
#ifndef __NOISE_ESTIMATE_H__
#define __NOISE_ESTIMATE_H__

namespace minc
{
  double noise_estimate(const minc::simple_volume<float>& input,double &mean_signal,bool gaussian=false,bool verbose=false);

};

#endif //__NOISE_ESTIMATE_H__