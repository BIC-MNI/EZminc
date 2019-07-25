/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2009 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include "pca_utils.h"

#include <fstream>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <cmath>
#include <gsl_glue.h>
#include <gsl/gsl_eigen.h>

#include <libgen.h> //for dirname

namespace minc
{
  
  void calculate_covariance(string_table& tbl, 
                            const volumes &mean, 
                            gsl_double_matrix& cov, 
                            int lo, int hi,
                            int skip,
                            bool normalize,
                            double dv
                            )
  {
    cov.resize(tbl.size()-skip,tbl.size()-skip);
    cov=0.0;
    minc::simple_volume<float> v1,v2;
    for(int k=lo;k<hi;k++)
    {
      for(int i=skip;i<tbl.size();i++)
      {
        std::cout<<i<<"\t"<<std::flush;
        minc_1_reader rdr1;
        rdr1.open(tbl[i][k].c_str());
        load_simple_volume<float>(rdr1,v1);
        if(mean.size()>=hi)
          v1-=mean[k];
        
        for(int j=skip;j<=i;j++)
        {
          if(i==j)
            v2=v1;
          else
          {
            minc_1_reader rdr2;
            rdr2.open(tbl[j][k].c_str());
            load_simple_volume<float>(rdr2,v2);
            
            if(mean.size()>=hi)
              v2-=mean[k];
          }
          v2*=v1;
          double ss=sum(v2);
          if(normalize) ss/=v2.size().vol();
          ss/=tbl.size()-skip;
          cov.set(i-skip,j-skip,cov.get(i-skip,j-skip)+ss*dv); //sum up all modalities
        }
      }
      std::cout<<std::endl;
    }
    //fill the upper triangle
    for(int i=0;i<(tbl.size()-skip);i++)
      for(int j=0;j<i;j++)
    {
      cov.set(j,i,cov.get(i,j));
    }
  }

  void calculate_covariance(string_table& tbl, 
                            const volumes &mean, 
                            const minc_byte_volume & mask,
                            gsl_double_matrix& cov, 
                            int lo, int hi,
                            int skip,
                            bool normalize,
                            double dv
                            )
  {
    cov.resize(tbl.size()-skip,tbl.size()-skip);
    cov=0.0;
    minc::simple_volume<float> v1,v2;
    for(int k=lo;k<hi;k++)
    {
      for(int i=skip;i<tbl.size();i++)
      {
        std::cout<<i<<"\t"<<std::flush;
        minc_1_reader rdr1;
        rdr1.open(tbl[i][k].c_str());
        load_simple_volume<float>(rdr1,v1);
        if(mean.size()>=hi)
          v1-=mean[k];

        for(int j=skip;j<=i;j++)
        {
          if(i==j)
            v2=v1;
          else
          {
            minc_1_reader rdr2;
            rdr2.open(tbl[j][k].c_str());
            load_simple_volume<float>(rdr2,v2);

            if(mean.size()>=hi)
              v2-=mean[k];
          }
          //v2*=v1;
          masked_mul(v2,v1,mask);
          double ss=sum(v2);

          if(normalize) 
            ss/=sum(mask);

          ss/=tbl.size()-skip;
          cov.set(i-skip,j-skip,cov.get(i-skip,j-skip)+ss*dv); //sum up all modalities
        }
      }
      std::cout<<std::endl;
    }
    //fill the upper triangle
    for(int i=0;i<(tbl.size()-skip);i++)
      for(int j=0;j<i;j++)
    {
      cov.set(j,i,cov.get(i,j));
    }
  }


  void calculate_covariance(string_table& tbl, 
                            const grids &mean, 
                            const minc_byte_volume & mask,
                            gsl_double_matrix& cov, 
                            int lo, int hi,
                            int skip,double dv
                           )
  {
    cov.resize(tbl.size()-skip,tbl.size()-skip);
    cov=0.0;
    minc_grid_volume v1,v2;
    for(int k=lo;k<hi;k++)
    {
      for(int i=skip;i<tbl.size();i++)
      {
        std::cout<<i<<"\t"<<std::flush;
        minc_1_reader rdr1;
        rdr1.open(tbl[i][k].c_str());
        load_simple_volume(rdr1,v1);
        if(mean.size()>=hi)
          v1-=mean[k];
        
        for(int j=skip;j<=i;j++)
        {
          if(i==j)
            v2=v1;
          else
          {
            minc_1_reader rdr2;
            rdr2.open(tbl[j][k].c_str());
            load_simple_volume(rdr2,v2);
            
            if(mean.size()>=hi)
              v2-=mean[k];
          }
          masked_mul(v2,v1,mask);
          double ss=sum(v2);
          //ss/=sum(mask);
          cov.set(i-skip,j-skip,cov.get(i-skip,j-skip)+ss*dv); //sum up all modalities
        }
      }
      std::cout<<std::endl;
    }
    //fill the upper triangle
    for(int i=0;i<(tbl.size()-skip);i++)
      for(int j=0;j<i;j++)
    {
      cov.set(j,i,cov.get(i,j));
    }
  }
  
  
  void calculate_covariance_cached(string_table& tbl, 
                                    const grids &mean, 
                                    const minc_byte_volume & mask,
                                    gsl_double_matrix& cov, 
                                    int lo, int hi,
                                    int skip,double dv
                                  )
  {
    
    std::vector<grids> inputs;
    inputs.resize(tbl.size());
    
    for(int j=skip;j<tbl.size();j++)
    {
      inputs.resize(tbl[0].size());
      for(int i=lo;i<hi;i++)
      {
        minc_1_reader rdr;
        rdr.open(tbl[j][i].c_str());
        load_simple_volume(rdr,inputs[j][i]);
        if(mean.size()>=hi)
          inputs[j][i]-=mean[i];
      }
    }
    
    
    cov.resize(tbl.size()-skip,tbl.size()-skip);
    cov=0.0;
    minc_grid_volume v2;
    for(int k=lo;k<hi;k++)
    {
      for(int i=skip;i<tbl.size();i++)
      {
        std::cout<<i<<"\t"<<std::flush;
        for(int j=skip;j<=i;j++)
        {
          v2=inputs[j][k];
          masked_mul(v2,inputs[i][k],mask);
          double ss=sum(v2);
          //ss/=sum(mask);
          cov.set(i-skip,j-skip,cov.get(i-skip,j-skip)+ss*dv); //sum up all modalities
        }
      }
      std::cout<<std::endl;
    }
    //fill the upper triangle
    for(int i=0;i<(tbl.size()-skip);i++)
      for(int j=0;j<i;j++)
    {
      cov.set(j,i,cov.get(i,j));
    }
  }
  
  void calculate_covariance(std::vector<gsl_double_matrix>& ipc,
                            const std::vector<int>& select,
                            const std::vector<double>& weights,
                            gsl_double_matrix& cov)
  {
    cov.resize(ipc[0].rows(), ipc[0].rows() );
    cov=0.0;
    
    for(int k=0;k<select.size();k++) {
      for(int i=0;i<ipc[0].rows();i++) {
        for(int j=0;j<=i;j++) 
        {
          double sum=0.0;
          for(int s=0;s<select[k];s++)
            sum+=ipc[k].get(s,i)*ipc[k].get(s,j);
          cov.set(i,j,cov.get(i,j)+sum*weights[k]/ipc[0].rows()); //sum up all modalities
        }
      }
    }
    
    //fill the upper triangle
    for(int i=0;i<ipc[0].rows();i++)
      for(int j=0;j<i;j++) 
    {
      cov.set(j,i,cov.get(i,j));
    }
  }
  
  int calculate_pca(gsl_double_matrix& cov,gsl_double_matrix& pc,gsl_double_vector &v,double threshold)
  {
    gsl_double_vector eigenvalue(cov.size(0));
    gsl_double_matrix eigvector(cov.size(0),cov.size(0));
        
    gsl_eigen_symmv_workspace * w=gsl_eigen_symmv_alloc(cov.size(0));
    gsl_eigen_symmv(cov, eigenvalue, eigvector, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eigenvalue, eigvector,GSL_EIGEN_SORT_ABS_DESC);
  
    double sum_eigenvalue=0.0;
    for(int i=0; i<cov.size(0); i++)    //select eigenvalue above the treshhold
    {
      sum_eigenvalue+=fabs(eigenvalue[i]);
    }
    
    double sub_sum=0.0;
    int n_select=0;
    
    for(int i=0; i< cov.size(0); i++)
    {
      sub_sum+=fabs(eigenvalue[i]);
      double eigen_ratio=sub_sum/sum_eigenvalue;
      if(eigen_ratio>threshold)
        break;
      n_select++;
    }
    
    pc=eigvector;
    v=eigenvalue;
    
    return n_select;
  }
}; //minc
