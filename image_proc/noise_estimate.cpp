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

static double noise_correct(double SIG,double NSIG)
{
  double SNR1=SIG/NSIG;

  for(int i=0;i<500;i++)
  {
    double SNR2 = sqrt(epsi(SNR1)*(1.0 + SIG*SIG  / (NSIG*NSIG) )-2.0);
    
    if( fabs(SNR1-SNR2) < 1e-9)
            break;
    SNR1 = SNR2;
  }    

  return sqrt((NSIG*NSIG / epsi(SNR1)));
}

double minc::noise_estimate(const minc::simple_volume<float>& input,double &mean_signal,bool gaussian,bool verbose)
{
  int debug=0;
  int maxiter=10;
  int hist_bins=2000;

	std::vector<simple_volume<float> > dwt;
	
	//Wavelet transform
	if(verbose) std::cout<<"Wavelet decomposition..."<<std::endl;
	dwt_forward(input,dwt);
	
	//detect object
	histogram<double> LLL_hist(hist_bins);
	build_histogram(LLL_hist,dwt[0]);
	
	minc_byte_volume LLL_mask(dwt[0].size());
	std::vector<double> LLL_mu;
	simple_k_means(LLL_hist,LLL_mu,2,10);
	
	double LLL_threshold=(LLL_mu[0]+LLL_mu[1])/2;
	
	if(verbose) std::cout<<"LLL threshold="<<LLL_threshold<<std::endl;
	
	for(int i=0;i<LLL_mask.c_buf_size();i++)
		LLL_mask.c_buf()[i]=dwt[0].c_buf()[i]>LLL_threshold ;
	
	
	double bkgr_mean=0.0;
	double bkgr_std=0.0;
	int bkgr_cnt=0;
	
	for(int z=0;z<input.dim(2);z++)
		for(int y=0;y<input.dim(1);y++)
			for(int x=0;x<input.dim(1);x++)
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
			
	// remove edges
	simple_volume<float> LLL_gmag(dwt[0].size());
	calc_gradient_mag(dwt[0],LLL_gmag,1,1,1);//TODO: compensate for nonuniform step size?
	
	histogram<double> gmag_hist(hist_bins);
	build_histogram(gmag_hist,LLL_gmag,LLL_mask);
	
	double gmag_median=gmag_hist.find_percentile( 0.5) ;
	
	if(verbose) std::cout<<"Median gradient magnitude="<<gmag_median<<std::endl;
	
	
	for(int i=0;i<LLL_mask.c_buf_size();i++)
		LLL_mask.c_buf()[i]=dwt[0].c_buf()[i]>LLL_threshold && LLL_gmag.c_buf()[i]<gmag_median;
	
	simple_volume<float> HHH_abs(dwt[0].size());
	
	for(int i=0;i<HHH_abs.c_buf_size();i++)
		HHH_abs.c_buf()[i]=fabs(dwt[7].c_buf()[i]);
		
	histogram<double> abs_HHH_hist(hist_bins);
	build_histogram(abs_HHH_hist,HHH_abs,LLL_mask);
	
	
	double nsig=abs_HHH_hist.find_percentile(0.5)/0.6745; //MAD estimator
	
	if(verbose) std::cout<<"Noise="<<nsig<<std::endl;
	
	/*
	if(verbose) std::cout<<"Blurring..."<<std::endl;
	simple_volume<float> blur(input.size());
	blur_volume(input,blur,false,false,false,3,3,3); //for now replace flat filter with gaussian
	
	//detect object
	histogram<double> blur_hist(hist_bins);
	build_histogram(blur_hist,blur);
	
	std::vector<double> blur_mu;
	simple_k_means(blur_hist,blur_mu,2,10);
	
	double blur_threshold=(blur_mu[0]+blur_mu[1])/2;
	if(verbose) std::cout<<"Signal threshold="<<    blur_threshold<<std::endl;      
	
	double mean_signal=0.0;
	int    cnt_signal=0;
	for(int i=0;i<input.c_buf_size();i++)
	{
		if(input.c_buf()[i]>blur_threshold)
		{
			mean_signal+=input.c_buf()[i];
			cnt_signal++;
		}
	}*/
	
	mean_signal=0.0;
	int    cnt_signal=0;
	
	for(int z=0;z<input.dim(2);z++)
		for(int y=0;y<input.dim(1);y++)
			for(int x=0;x<input.dim(1);x++)
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
  
  //if(verbose)
  //  std::cout<<"Noise="<<nsig<<std::endl;
		
	return nsig_corr;
}
