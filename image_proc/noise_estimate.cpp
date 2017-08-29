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
#include <iostream>
#include <getopt.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include "minc_histograms.h"
//#include "pca_utils.h"
#include "fftw_blur.h"
#include "dwt_utils.h"
#include "noise_estimate.h"
#include <gsl/gsl_sf_bessel.h>

using namespace minc;

static double epsi(double SNR)
{
  if (SNR > 37) return 1.0;
  
  return 2.0 + SNR*SNR - M_PI/8 * exp(-(SNR*SNR/2.0))*pow((2.0+SNR*SNR)*gsl_sf_bessel_I0(SNR*SNR/4.0) + SNR*SNR*gsl_sf_bessel_I1(SNR*SNR/4.0),2);

}

static double noise_correct(double sig,double nsig)
{
  double snr1=sig/nsig;

  for(int i=0;i<500;i++)
  {
    double snr2 = sqrt(epsi(snr1)*(1.0 + sig*sig  / (nsig*nsig) )-2.0);
    
    if( fabs(snr1-snr2) < 1e-9)
            break;
    snr1 = snr2;
  }    

  return sqrt((nsig*nsig / epsi(snr1)));
}

double minc::noise_estimate(const minc::simple_volume<float>& input,double &mean_signal,bool gaussian,bool verbose,int hist_bins,const minc::minc_byte_volume& mask)
{
  int debug=0;

	std::vector<simple_volume<float> > dwt;
	
	//Wavelet transform
	if(verbose) std::cout<<"Wavelet decomposition..."<<std::endl;
	dwt_forward(input,dwt);
	
	//detect object
	histogram<double> LLL_hist(hist_bins);
	build_histogram(LLL_hist,dwt[0]);
	
	minc_byte_volume LLL_mask(dwt[0].size());
  
  if(mask.empty()) {
    std::vector<double> LLL_mu;
    simple_k_means(LLL_hist,LLL_mu,2,10);
    
    float LLL_threshold=(float)(LLL_mu[0]+LLL_mu[1])/2;
    
    if(verbose) std::cout<<"LLL threshold="<<LLL_threshold<<std::endl;
    
    for(size_t i=0;i<LLL_mask.c_buf_size();i++)
      LLL_mask.c_buf()[i]=(dwt[0].c_buf()[i]>LLL_threshold?1:0);
  } else {
    /*downsample mask */
    LLL_mask=0;
    for(size_t z=0;z<input.dim(2)/2;z++)
      for(size_t y=0;y<input.dim(1)/2;y++)
        for(size_t x=0;x<input.dim(0)/2;x++)
    {
      LLL_mask.set(x,y,z,mask.get(x*2,y*2,z*2));
    }
  }
	
	double bkgr_mean=0.0;
	double bkgr_std=0.0;
	size_t bkgr_cnt=0;
	
	for(size_t z=0;z<input.dim(2);z++)
		for(size_t y=0;y<input.dim(1);y++)
			for(size_t x=0;x<input.dim(0);x++)
	{
		if(!LLL_mask.get(x/2,y/2,z/2))
		{
			double v=input.get(x,y,z);
			bkgr_mean+=v;
			bkgr_std+=v*v;
			bkgr_cnt++;
		}
	}
	
	if(bkgr_cnt>0)
	{
		bkgr_mean/=bkgr_cnt;
		bkgr_std/=bkgr_cnt;
		bkgr_std-=bkgr_mean*bkgr_mean;
		bkgr_std=sqrt(bkgr_std);
		if(verbose) {
			std::cout<<"Background mean="<<bkgr_mean<< " std="<<bkgr_std<<std::endl;
		}
	}
			
	// remove edges (areas with high magnitude of the gradients)
	simple_volume<float> LLL_gmag(dwt[0].size());
	calc_gradient_mag(dwt[0],LLL_gmag,1,1,1);//TODO: compensate for nonuniform step size?
	
	histogram<double> gmag_hist(hist_bins);
	build_histogram(gmag_hist,LLL_gmag,LLL_mask);
	
	double gmag_median=gmag_hist.find_percentile( 0.5) ;
	
	if(verbose) std::cout<<"Median gradient magnitude="<<gmag_median<<std::endl;
	
  //exclude areas with high gradients from mask
	for(int i=0;i<LLL_mask.c_buf_size();i++)
		LLL_mask.c_buf()[i]=LLL_mask.c_buf()[i] && LLL_gmag.c_buf()[i]<gmag_median;
	
	simple_volume<float> HHH_abs(dwt[0].size());
	
	for(int i=0;i<HHH_abs.c_buf_size();i++)
		HHH_abs.c_buf()[i]=fabs(dwt[7].c_buf()[i]);
		
	histogram<double> abs_HHH_hist(hist_bins);
	build_histogram(abs_HHH_hist,HHH_abs,LLL_mask);
	
	
	double nsig=abs_HHH_hist.find_percentile(0.5)/0.6745; //MAD estimator
	
	if(verbose) std::cout<<"Noise="<<nsig<<std::endl;
	
	mean_signal=0.0;
	size_t    cnt_signal=0;
	
	for(size_t z=0;z<input.dim(2);z++)
		for(size_t y=0;y<input.dim(1);y++)
			for(size_t x=0;x<input.dim(0);x++)
	{
		if(LLL_mask.get(x/2,y/2,z/2))
		{
			mean_signal+=input.get(x,y,z);
			cnt_signal++;
		}
	}
	
	if(cnt_signal) mean_signal/=cnt_signal;
  
	if(verbose) std::cout<<"Signal="<<mean_signal<<std::endl;
  double nsig_corr=nsig;
  
  if(!gaussian)
  {
    if(verbose) std::cout<<"Correcting SNR based on Koay method"<<std::endl;
    nsig_corr=noise_correct(mean_signal,nsig);
  }
  
  return nsig_corr;
}

// kate: indent-mode cstyle; indent-width 2; replace-tabs on; tab-width 2
