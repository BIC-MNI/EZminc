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
#ifndef __GSL_NLFIT_H__
#define __GSL_NLFIT_H__

#include "gsl_glue.h"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

namespace minc
{

  class gsl_cost_function
  {
  public:
    
    //! override this
    virtual int expb_f(const gsl_double_vector& x, gsl_double_vector& f)
    {
      return GSL_FAILURE;
    }
    
    virtual int expb_df(const gsl_double_vector& x, gsl_double_matrix& J)
    {
      return GSL_FAILURE;
    }
    
    virtual int expb_fdf(const gsl_double_vector& x, gsl_double_vector& f, gsl_double_matrix& J)
    {
     if(expb_f(x,f)==GSL_SUCCESS && expb_df(x,J)==GSL_SUCCESS)
       return GSL_SUCCESS;
     return GSL_FAILURE;
    }
    
    virtual ~gsl_cost_function()
    {
    }
  };
  
  class gsl_nlfit
  {
  protected:
    gsl_multifit_fdfsolver   *_solver;
    gsl_multifit_function_fdf _fun_f;
    int _status;
    void clean()
    {
      if(_solver)
        gsl_multifit_fdfsolver_free(_solver);
      _solver=NULL;
    }
    static int expb_f   (const gsl_vector * x, void *data, gsl_vector * f);
    static int expb_df  (const gsl_vector * x, void *data, gsl_matrix * J);
    static int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);
    double _x_delta,_dx_delta;
    
    mutable gsl_double_vector _current;
    gsl_double_vector _grad;
    mutable gsl_double_matrix _cur_j;
    gsl_double_matrix _cov;
    gsl_double_vector _error;
  public:
    
    const gsl_double_vector& get_error(void) const
    {
      return _error;
    }
    
    gsl_nlfit():_solver(NULL),_x_delta(1e-4),_dx_delta(1e-4)
    {
    }
    
    ~gsl_nlfit()
    {
      clean();
    }
    
    void set_tol(double d_x,double d_dx)
    {
      _x_delta=d_x;
      _dx_delta=d_dx;
    }
    
    void init(int n, int p, gsl_cost_function* fun);
    
    void start(gsl_double_vector& x)
    {
      gsl_multifit_fdfsolver_set(_solver,&_fun_f,x);
    }
    
    bool iterate(void)
    {
      _status = gsl_multifit_fdfsolver_iterate (_solver);
      if(_status) return false;
      gsl_multifit_gradient(_solver->J,_solver->f,_grad);
      
      _status = gsl_multifit_test_delta(_solver->dx, _solver->x, _dx_delta, _x_delta);
      //_status= gsl_multifit_test_gradient(_grad,1e-6);
      return _status==GSL_CONTINUE;
    }

    const gsl_double_vector& solution(void) const
    {
      //_current.attach(_solver->x);
      return _current;
    }

    const gsl_double_matrix& J(void) const
    {
      //_cur_j.attach(_solver->J);
      return _cur_j;
    }

    const gsl_double_matrix& calc_covariance(void);

    const gsl_double_vector& calc_error(void);

    double rms(void) const;
    
  };
};

#endif //__GSL_NLFIT_H__
