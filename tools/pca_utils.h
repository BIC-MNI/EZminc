#ifndef __PCA_UTILS_H__
#define __PCA_UTILS_H__

#include "utils.h"
#include <gsl_glue.h>

namespace minc
{
  //! calculate covariance matrix, only  pair of files is loaded into memory
  void calculate_covariance(string_table& tbl,const volumes &mean,
                            gsl_double_matrix& cov,int lo,int hi,
                            int skip=0,bool normalize=false, double dv=1.0);
                            
  void calculate_covariance(string_table& tbl, const volumes &mean, 
                            const minc_byte_volume & mask,
                            gsl_double_matrix& cov,int lo, int hi,
                            int skip=0,bool normalize=false, double dv=1.0);
 
  
  //! calculate covariance matrix, only  pair of files is loaded into memory
  void calculate_covariance(string_table& tbl,const grids &mean,
                            const minc_byte_volume& mask,
                            gsl_double_matrix& cov,int lo,int hi,
                            int skip=0, double dv=1.0);
  
  //! calculate covariance matrix, all files are loaded into memory
  void calculate_covariance_cached(string_table& tbl,const grids &mean,
                            const minc_byte_volume& mask,
                            gsl_double_matrix& cov,int lo,int hi,
                            int skip=0, double dv=1.0);
  
  //! calculate covariance matrix, from the multimodal principal components, using only select components
  void calculate_covariance(std::vector<gsl_double_matrix>& ipc,
                            const std::vector<int>& select,
                            const std::vector<double>& weights,
                            gsl_double_matrix& cov);

  //! calculate PCA from a covariance matrix
  int calculate_pca(gsl_double_matrix& cov,gsl_double_matrix& pc,gsl_double_vector &v,double threshold);
  
  //! calculate number of PCs to represent a given fraaction of variance (1 - 100%)
  int calc_selected(string_table& tbl,double threshold);
  
};
#endif //__PCA_UTILS_H__
