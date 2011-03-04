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


#include "fftw_blur.h"
#include <fftw3.h>
#include <math.h>
#include <iostream>
namespace minc
{
  

  void calculate_gaussian(std::vector<complex_f>& res,double sigma,bool dx)
  {
    double pi=2*atan2(1.0,0.0);
  
    double a=1.0/(2.0*sigma*sigma);
  
    int n=res.size();
  
    double sum=0.0;
    for(int i=0;i<n;i++)
    {
      if(i<(n/4)||i>=(3.0*n/4.0))
      {
        double k=(i<n/2.0?i:i-n);
        
        sum+=exp(-k*k/(2.0*sigma*sigma));

        if(dx)
          res[i]=-k/(sigma*sigma)*exp(-k*k/(2.0*sigma*sigma));
        else
          res[i]=exp(-k*k/(2.0*sigma*sigma));
      } else {
        res[i]=0.0;
      }
    }
  
  //if(!dx)
    for(int i=0;i<n;i++)
    {
      res[i]/=sum;
    }
  }
  

  void blur_volume(simple_volume<float> &in,simple_volume<float> &out, bool dx,bool dy, bool dz, double fx,double fy,double fz)
  {
    double sx=fabs(fx/(2.0*sqrt(2*log(2.0))));
    double sy=fabs(fy/(2.0*sqrt(2*log(2.0))));
    double sz=fabs(fz/(2.0*sqrt(2*log(2.0))));
    //std::cout<<sx<<" "<<sy<<" "<<sz<<std::endl;
    fftwf_plan px,ipx,py,ipy,pz,ipz;
    std::vector<complex_f> _fx(in.dim(0)*2),_gx(in.dim(0)*2);
    std::vector<complex_f> _fy(in.dim(1)*2),_gy(in.dim(1)*2);
    std::vector<complex_f> _fz(in.dim(2)*2),_gz(in.dim(2)*2);
    
    px=fftwf_plan_dft_1d(_fx.size(), 
                        reinterpret_cast<fftwf_complex*>(&_fx[0]), 
                        reinterpret_cast<fftwf_complex*>(&_fx[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    ipx=fftwf_plan_dft_1d(_fx.size(),
                          reinterpret_cast<fftwf_complex*>(&_fx[0]), 
                          reinterpret_cast<fftwf_complex*>(&_fx[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
  
    py=fftwf_plan_dft_1d(_fy.size(), 
                        reinterpret_cast<fftwf_complex*>(&_fy[0]), 
                        reinterpret_cast<fftwf_complex*>(&_fy[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    ipy=fftwf_plan_dft_1d(_fy.size(),
                          reinterpret_cast<fftwf_complex*>(&_fy[0]), 
                          reinterpret_cast<fftwf_complex*>(&_fy[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
  
    pz=fftwf_plan_dft_1d(_fz.size(), 
                        reinterpret_cast<fftwf_complex*>(&_fz[0]), 
                        reinterpret_cast<fftwf_complex*>(&_fz[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    ipz=fftwf_plan_dft_1d(_fz.size(),
                        reinterpret_cast<fftwf_complex*>(&_fz[0]), 
                        reinterpret_cast<fftwf_complex*>(&_fz[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
    
    //precalculate gaussians in frequency domain 
    calculate_gaussian(_fx,sx,dx);
    calculate_gaussian(_fy,sy,dy);
    calculate_gaussian(_fz,sz,dz);
    fftwf_execute(px);
    fftwf_execute(py);
    fftwf_execute(pz);
    _gx=_fx;
    _gy=_fy;
    _gz=_fz;
    
    int i,j,k;
    
    //perform calculation in Z
    if(sz>0.0)
    {
      for(i=0;i<in.dim(0);i++)
        for(j=0;j<in.dim(1);j++)
      {
        for(k=0;k<in.dim(2);k++)
        {
          _fz[k]=in(i,j,k);
          _fz[k+in.dim(2)]=0.0;
        }
        fftwf_execute(pz);
          
        for(k=0;k<_fz.size();k++)
          _fz[k]*=_gz[k];
            
        fftwf_execute(ipz);
        for(k=0;k<in.dim(2);k++)
        {
          out(i,j,k)=_fz[k].real()/(in.dim(2)*2);
        }
      }
    } else  {
      for(i=0;i<in.dim(0);i++)
        for(j=0;j<in.dim(1);j++)
          for(k=0;k<in.dim(2);k++)
            out(i,j,k)=in(i,j,k); 
    }
        
    if(sy>0.0)
    {
      for(i=0;i<in.dim(0);i++)
        for(k=0;k<in.dim(2);k++)
      {
        for(j=0;j<in.dim(1);j++)
        {
          _fy[j]=out(i,j,k);
          _fy[j+in.dim(1)]=0.0;
        }
        fftwf_execute(py);
          
        for(j=0;j<_fy.size();j++)
          _fy[j]*=_gy[j];
            
        fftwf_execute(ipy);
        for(j=0;j<in.dim(1);j++)
        {
          out(i,j,k)=_fy[j].real()/(in.dim(1)*2);
        }
      }
    } 
    
    if(sx>0.0)
    {
      for(j=0;j<in.dim(1);j++)
        for(k=0;k<in.dim(2);k++)
      {
        for(i=0;i<in.dim(0);i++)
        {
          _fx[i]=out(i,j,k);
          _fx[i+in.dim(0)]=0.0;
        }
        fftwf_execute(px);
          
        for(i=0;i<_fx.size();i++)
          _fx[i]*=_gx[i];
            
        fftwf_execute(ipx);
        for(i=0;i<in.dim(0);i++)
        {
          out(i,j,k)=_fx[i].real()/(in.dim(0)*2);
        }
      }
    }
      
    fftwf_destroy_plan(px);fftwf_destroy_plan(ipx);
    fftwf_destroy_plan(py);fftwf_destroy_plan(ipy);
    fftwf_destroy_plan(pz);fftwf_destroy_plan(ipz);
  }
  
  void calc_gradient(simple_volume<float> &in,minc_grid_volume &out, double fx,double fy,double fz)
  {
    simple_volume<float> vdx(in); 
    simple_volume<float> vdy(in); 
    simple_volume<float> vdz(in); 
    
    blur_volume(in,vdx,1,0,0,fx,fy,fz);
    blur_volume(in,vdy,0,1,0,fx,fy,fz);
    blur_volume(in,vdz,0,0,1,fx,fy,fz);
    
    out.resize(in.size());
    
    for(int i=0;i<out.c_buf_size();i++) //calculating magnitude
      out.c_buf()[i]=IDX<float>(vdx.c_buf()[i],vdy.c_buf()[i],vdz.c_buf()[i]);
    
  }
  
  void calc_gradient_mag(simple_volume<float> &in,simple_volume<float> &out, double fx,double fy,double fz)
  {
    
    minc_grid_volume tmp(in.size());
    calc_gradient(in,tmp,fx,fy,fz);
    
    for(int i=0;i<out.c_buf_size();i++) //calculating magnitude
    {
      out.c_buf()[i]=sqrt(tmp.c_buf()[i].mod2());
    }
  }
  
};
