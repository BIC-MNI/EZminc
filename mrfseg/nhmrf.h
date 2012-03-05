// ******************************************************
// Functions for SVPA based non homogeneous MRFs 
// MRFSEG + GAMIXTURE software bundle: Voxel classification in volumetric 
// images using genetic algorithms and Markov random fields.
// The genetic algorithm is described in: 
// [1]  J. Tohka , E. Krestyannikov, I.D. Dinov , A. MacKenzie-Graham,
// D.W. Shattuck , U. Ruotsalainen, and A.W. Toga. 
// Genetic algorithms for finite mixture model based voxel classification
// in neuroimaging.  
// IEEE Transactions on Medical Imaging  , 26(5):696 - 711, 2007.
// (C) Jussi Tohka, 2005 - 2007 
// Laboratory of NeuroImaging, Department of Neurology, UCLA Medical School, CA, USA
// Institute of Signal Processing, Tampere University of Technology, Finland  
// Version 1.1
// e-mail jussi.tohka@tut.fi
// *****************************************************
// ******************************************************
// Permission to use, copy, modify, and distribute this software 
// for any purpose and without fee is hereby
// granted, provided that the above copyright notice appears in all
// copies.  The author, Tampere University of Technology and University of
// California, Los Angeles make no representations
// about the suitability of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
// *****************************************************

#ifndef NHMRF_H
#define NHMRF_H

#include <vector>

#include "generic_image.h"
#include "atlasspec.h"

#define VERY_SMALL 1E-16
#define INTERVALS 100

#define DEFAULT_PVESTART 0.0
#define DEFAULT_PVEEND 1.0


class MixtureSpec { 
  public: //for now
  AtlasSpec* patlas;              // remember that forbidden labels will mean that the corresponding probability is zero
  std::vector<float> prob;        //  prior probabilities
  std::vector<float> mu;          //  means, 
  std::vector<float> sigma2;      //  variances
  
  void allocateMixtureSpec(AtlasSpec& atlas)
  {
    patlas = &atlas;
    prob.resize((atlas.n) * (atlas.numberOfLabels));
    mu.resize((atlas.n) * (atlas.numberOfLabels));
    sigma2.resize((atlas.n) * (atlas.numberOfLabels));   
  }
  
  void copyMixtureSpec(const MixtureSpec& a)
  {
    prob=a.prob;
    mu=a.mu;
    sigma2=a.sigma2;
  }
  
  float getProb(int region,int label) const
  {
    return(prob[(patlas->numberOfLabels)*region + label]);  
  };
  
  void putProb(int region,int label, float val)
  {
    prob[(patlas->numberOfLabels)*region + label] = val;  
  };
  
  float getMu(int region,int label)
  {
    return(mu[(patlas->numberOfLabels)*region + label]);  
  };
  
  void putMu(int region,int label, float val)
  {
    mu[(patlas->numberOfLabels)*region + label] = val;  
  };
  
  float getSigma2(int region,int label)
  {
    return(sigma2[(patlas->numberOfLabels)*region + label]);  
  };

  void putSigma2(int region,int label, float val)
  {
    sigma2[(patlas->numberOfLabels)*region + label] = val;  
  };

  void printMixture();
};

// Normalizes n probalities to sum to one
inline void normalize(float* pval, int n)

{ float sca = 0.0;
  int i;

  for(i = 0;i < n;i++) {
    sca = sca + pval[i];
  }
  if(fabs(sca) >  VERY_SMALL) {       // To avoid divisions by zero 
    for(i = 0;i < n;i++) {
      pval[i] = pval[i]/sca;
    }
  }
}

// Finds maximum argument out of the n possibilities
// moved to analyze.h
// inline char maxArg(float* pval,int n)
// {
//  float maximum;
//  int i,index;
//  
//  maximum = pval[0];
//  index = 0;
//  for(i = 1;i < n;i++) {
//    if(pval[i] > maximum) {
//      index = i;
//      maximum = pval[i];
//    }
//  }
//  return( (char) index);
// }

inline void collectValuesFromImagePP(const std::vector<FloatImage>& imagePP,float* collectHere,
                                     int x, int y, int z,int n) 
{
  int i;
  for(i = 0;i < n;i++) {
    collectHere[i] =imagePP[i].getVoxelValue(x,y,z);
  } 
}

// Computes likelihood of value given parameters mean and variance. 
// Returns the likelihood.

inline float computeGaussianLikelihood(float value, float mean , float var)

{ 
  return(exp(-((value - mean) *(value - mean))/(2 * var))/(sqrt(2 * M_PI * var)));

}

// Computes the likelihoods for the mixed classes. Returns the likelihood.
// var1,var2 are the variances of pdfs representing pure classes.
// So the model for the variable y (representing the 
// intensity value) that is composed of t * tissue1 and (1 - t)* tissue2 becomes :

// y = t*x1 + (1 - t)*x2 ,
// x1 ~ N(mean1,var1) , x2 ~ N(mean2,var2) 

inline float computeMarginalizedLikelihood(float value, float mean1 , float mean2, 
                                            float var1, float var2, 
					      unsigned int nof_intervals, float start, float end)

{ 
  float lh, tmean , tvar, t, interval_len;
  int i;  
  
  interval_len = (float) (end - start) / nof_intervals;
  lh = 0;
  for(i = 0; i < nof_intervals; i++) {
    t = (i + 0.5) * interval_len + start;
    tmean = t * mean1 + ( 1 - t ) * mean2;
    tvar = pow(t,2) * var1 + pow((1 - t),2) * var2;
    lh = lh + computeGaussianLikelihood(value, tmean,tvar) / nof_intervals;
  }
  return(lh);
}

inline float secondOrderGibbs(char testLabel,const LabelImage& labels,float* mrfConstants,
                              int x, int y,int z,float* distanceLookup, float beta)
{
  int i,j,k;
  float exponent = 0.0;
  float tmp;

  for(i = (-1);i < 2;i++) {
    for(j = (-1);j < 2;j++) {
      for(k = (-1);k < 2;k++) {
        if(! ((i == 0) && (j == 0) && (k == 0)) ) {
          tmp = mrfConstants[labels.getSafeLabelValue(x + i,y + j, z + k)]/distanceLookup[ (i + 1) * 9  + (j + 1) * 3 + (k + 1) ];
          exponent = exponent + tmp;
	  
        }
      }
    }
  }  
  return(exp(- (beta * exponent)));
}


int readMixtureParameters(char* filename,AtlasSpec& atlas,MixtureSpec& mixture);

int writeMixtureParameters(char* filename,MixtureSpec& mixture,bool overwrite);

int computeVoxelLikelihood(MixtureSpec& mixture,const FloatImage& img,const FloatImage& mask,
    const std::vector<FloatImage>& atlasImages,std::vector<FloatImage>& labelLikelihoods);

int computeMRF(LabelImage& labels,MixtureSpec& mixture,const FloatImage& mask,
    const std::vector<FloatImage>& labelLikelihoods, const std::vector<FloatImage>& atlasImages, 
    float beta1, float beta2,int maxIterations, bool verbose);

int computeGibbs(LabelImage& labels,MixtureSpec& mixture, const FloatImage& mask,const std::vector<FloatImage>& labelLikelihoods, 
    const std::vector<FloatImage>& atlasImages, float beta1,float beta2,int maxIterations, bool verbose );

int convertPVElabels(LabelImage& crispLabels, const LabelImage& pveLabels, const FloatImage& img, 
    const std::vector<FloatImage>& atlasImages, MixtureSpec& mixture);


#endif
