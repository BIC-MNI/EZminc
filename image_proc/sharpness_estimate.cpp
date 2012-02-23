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
#include <iostream>
#include <getopt.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include "minc_histograms.h"
#include "fftw_blur.h"
#include "sharpness_estimate.h"

using namespace minc;


double minc::sharpness_estimate(const minc::simple_volume<float>& input,const minc::simple_volume<unsigned char>& mask)
{
  int debug=0;
  int hist_bins=2000;

			
	// remove edges
	simple_volume<float> gmag(input.size());
	calc_gradient_mag(input,gmag,1,1,1);//TODO: compensate for nonuniform step size?
	
	histogram<double> gmag_hist(hist_bins);
	build_histogram(gmag_hist,gmag,mask);
	
	double gmag_median=gmag_hist.find_percentile( 0.5) ;
	
	return gmag_median;
}
