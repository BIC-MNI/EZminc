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

double kl_distance(const histogram<double>& sample1,const histogram<double>& sample2)
{
  //let's include the whole range 
  double _min=std::min(sample1.min(),sample2.min());
  double _max=std::max(sample1.max(),sample2.max());
  double _range=_max-_min;
  
  //now we are going to iterate through samples
  double distance=0.0;
  int _size=std::max(sample1.size(),sample2.size());
  double _bin=_range/_size;
  double v=_min+_bin/2.0;
  
  for(int i=0;i<_size;i++,v+=_bin)
  {
    double _P=0.0;
    if(v>=sample1.min() && v<=sample1.max()) _P=sample1[v];
    double _Q=0.0;
    if(v>=sample2.min() && v<=sample2.max()) _Q=sample2[v];
    
    if(_P>0.0 && _Q>0.0)//TODO: epsilon?
      distance+=_P*log(_P/_Q);
  }
  return distance;
}

double ks_distance(const histogram<double>& sample1,const histogram<double>& sample2)
{
  //let's include the whole range 
  double _min=std::min(sample1.min(),sample2.min());
  double _max=std::max(sample1.max(),sample2.max());
  double _range=_max-_min;
  
  histogram<double> _s1(sample1);
  histogram<double> _s2(sample2);
  
  _s1.convert_to_commulative();
  _s2.convert_to_commulative();
  
  int _size=std::max(_s1.size(),_s2.size())*2;
  double _bin=_range/_size;
  
  double dist=0.0;
  double v=0.0;
  for(int i=0;i<_size;i++,v+=_bin)
  {
    double d=fabs(_s1[v]-_s2[v]);
    if(d>dist) dist=d;
  }
  return dist;
}


// Code adapted from "Numerical Recipes in C" 
static double probks(double alam)
{
  const double EPS1=0.001;
  const double EPS2=1e-8;
  
  double a2,fac=2.0,sum=0.0,term,termbf=0.0;
  a2 = -2.0*alam*alam;
  
  for (int j=1;j<=100;j++) {
    term=fac*exp(a2*j*j);
    sum += term;
    if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) 
      return sum;
    fac = -fac;
    termbf=fabs(term);
  }
  return 1.0;
}

double ks_significance(double dist, double n1,double n2)
{
  double en=sqrt(n1*n2/(n1+n2));
  return probks((en+0.12+0.11/en)*dist);
}

};