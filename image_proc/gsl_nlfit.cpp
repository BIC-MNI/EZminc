/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: interface to LSQ nonlinear optimization solver from GSL
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
#include "gsl_nlfit.h"

namespace minc
{
  void gsl_nlfit::init(int n, int p, gsl_cost_function* fun)
  {
    _solver = gsl_multifit_fdfsolver_alloc (gsl_multifit_fdfsolver_lmsder, n, p);
    _fun_f.f =   &expb_f;
    _fun_f.df  = &expb_df;
    _fun_f.fdf = &expb_fdf;
    _fun_f.n = n;
    _fun_f.p = p;
    _fun_f.params = fun;
      
    _current.attach(_solver->x);
    _cur_j.attach(_solver->J);
    _cov.resize(_fun_f.p,_fun_f.p);
      
    _grad.resize(p);
    _error.resize(p);
  }
  
  int gsl_nlfit::expb_f(const gsl_vector* x, void *data, gsl_vector * f)
  {
    gsl_double_vector _x((gsl_vector* )x);//this is a hack
    gsl_double_vector _f(f);//this is a hack
    return ((gsl_cost_function*)data)->expb_f(_x, _f);
  }
  
  int gsl_nlfit::expb_df(const gsl_vector * x, void *data, gsl_matrix * J)
  {
    gsl_double_vector _x((gsl_vector* )x);//this is a hack
    gsl_double_matrix _J(J);
    return ((gsl_cost_function*)data)->expb_df(_x,_J);
  }
  
  int gsl_nlfit::expb_fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
  {
    gsl_double_vector _x((gsl_vector* )x);//this is a hack
    gsl_double_vector _f(f);//this is a hack
    gsl_double_matrix _J(J);
    return ((gsl_cost_function*)data)->expb_fdf(_x,_f,_J);
  }
  
  const gsl_double_matrix& gsl_nlfit::calc_covariance(void)
  {
      //_cur_j.attach(_solver->J);
    gsl_multifit_covar(_solver->J, FLT_EPSILON, _cov);  
    return _cov;
  }
 
  const gsl_double_vector& gsl_nlfit::calc_error(void)
  {
    calc_covariance();
    double chi=gsl_blas_dnrm2(_solver->f);
    double c = chi / sqrt(_fun_f.n - _fun_f.p);
    
  
    for(int i=0;i<_fun_f.p;i++)
    {
      _error.set(i,sqrt(_cov.get(i,i))*c);
    }
    return _error;
  }
  
  double gsl_nlfit::rms(void) const
  {
    double sum=0.0;
    gsl_cost_function* cost=(gsl_cost_function*)_fun_f.params;
    gsl_double_vector f(_fun_f.n);
    cost->expb_f(_current,f);
    for(int i=0;i<_fun_f.n;i++)
    {
      sum+=f.get(i)*f.get(i);
    }
    return sqrt(sum/_fun_f.n);
  }
  
  
}; //minc
