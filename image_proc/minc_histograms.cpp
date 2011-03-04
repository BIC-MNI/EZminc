/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: image histogram manipulation routines
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
#include "minc_histograms.h"

namespace minc
{

//! estimate sample mu using discrete classes
void estimate_mu(const histogram<double>& input,
                 const std::vector<int>& cls,
                 std::vector<double>&  mu )
{
  int number_of_classes=mu.size();
  int i,j,k;
  std::vector<double>    _counts(number_of_classes,0);
  
  for(k=0;k<number_of_classes;k++)
    mu[k]=0;
    
  double _count=0;
  
  //1 calculate means
  for(j=0;j<cls.size();j++)
  {
    _count+=input[j];
    
    int _c=cls[j];
    if(!_c || _c>number_of_classes) continue; //only use classified voxel
    _c--;
    
    _counts[_c]+=input[j];
    mu[_c]+=input[j]*input.value(j);
  }
  
  if(!_count)
    REPORT_ERROR("No voxels defined in ROI!");
  
  for(k=0;k<number_of_classes;k++)
  {
    for(k=0;k<number_of_classes;k++)
    {
      if(_counts[k]>0)
        mu[k]/=_counts[k];
    }
  }
}

//!  apply k-means classify to histogram
void apply_k_means_classify(const histogram<double>& input,
                    const std::vector<double>&  mu,
                    std::vector<int>& cls )
{
  int number_of_classes=mu.size();
  int i,j,k;
  
  //calculate all the classes
  for(j=0;j<input.size();j++)
  {
    int best_k=0;
    double best_dist=0;
    
    for(k=0;k<number_of_classes;k++)
    {
      double dist=fabs(input.value(j)-mu[k]);
      if(dist<best_dist|| best_k==0)
      {
        best_dist=dist;
        best_k=k+1;
      }
    }
    cls[j]=best_k;
  }
}

//!  apply k-means classify to volume
void apply_hard_classify(simple_volume<float> & input,
                          minc_byte_volume    &mask,
                          std::vector<double>  mu,
                          minc_byte_volume& cls
                        )
{
  int number_of_classes=mu.size();
  int i,j,k;
  
  if(cls.size()!=mask.size())
    REPORT_ERROR("Mask is wrong size");
  
  //calculate all the classes
  for(j=0;j<input.c_buf_size();j++)
  {
    if(mask.c_buf()[j]) 
    {
      int best_k=0;
      double best_dist=0;
      for(k=0;k<number_of_classes;k++)
      {
        double dist=fabs(input.c_buf()[j]-mu[k]);
        if(dist<best_dist|| best_k==0)
        {
          best_dist=dist;
          best_k=k+1;
        }
      }
      cls.c_buf()[j]=best_k;
    } else {
      cls.c_buf()[j]=0;//unclassified?
    }
  }
}

void simple_k_means( histogram<double>& hist, std::vector<double>& mu,int k_means,int maxiter)
{
  std::vector<int>    cls(hist.size(),0);
  
  if(mu.size()!=k_means)
  {
  
    mu.resize(k_means);
    //initializing k-means using uniform distribution
    double vol_min=hist.find_percentile(0.001);
    double vol_max=hist.find_percentile(0.999);
        
    for(int j=0;j<k_means;j++)
        mu[j]= vol_min+(vol_max-vol_min)*j/(double)(k_means-1);
  }
        
  for(int iter=0;iter<maxiter;iter++)
  {
    apply_k_means_classify(hist,mu,cls);
    estimate_mu(hist,cls,mu);
  }
}

};