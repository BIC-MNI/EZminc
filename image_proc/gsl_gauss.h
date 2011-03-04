/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: interface to LSQ (Gauss) optimization solver from GSL
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
#ifndef __GSL_GAUSS_H__
#define __GSL_GAUSS_H__

#include "gsl_glue.h"
#include <gsl/gsl_linalg.h>

namespace minc 
{
	template < class T, class B > class MNK_Gauss
	{
	protected:
		//basis_function basis;
		int _basis_order;
		gsl_double_matrix _alpha;
		gsl_double_vector _beta;
		gsl_double_vector _x;
  
		gsl_double_matrix _V;
		gsl_double_vector _S;
		gsl_double_vector _work;
  
		B _basis;
  
    
	public:
		typedef std::vector < double > valuesVec;
		typedef std::vector < float >  valuesVecF;
    
		double fit (const valuesVec & coeff, T pnt) const
		{
			double val = 0.0;
			for (int i = 0; i < dim (); i++)
				val += _basis (i, pnt) * coeff[i];
			return val;
		}
    
		double fit (const valuesVec& basis,const valuesVec & coeff) const
		{
			double val = 0.0;
			for (int i = 0; i < dim (); i++)
				val += basis[i] * coeff[i];
			return val;
		}
	
		double fit (const valuesVec & coeff, T pnt, T min, T max) const
		{
			double val = fit (coeff, pnt);
			if (val > max)
				return max;
			if (val < min)
				return min;
			return val;
		}
	
		double basis (int i, T val)
		{
			return _basis (i, val);
		}
	
		int dim (void) const
		{
			return _basis_order;
		}
	
		//! basic construcor
		MNK_Gauss (int n,		//!< order of aproximation
							 B & basis	//!< basis functions
							):_basis_order (n), 
                _basis(basis), 
                _alpha (_basis_order, _basis_order),
                _beta (_basis_order),_V(_basis_order,_basis_order),
                _S(_basis_order),_work(_basis_order),_x(_basis_order)
		{
      clear();
		};
	
		//! basic construcor
		MNK_Gauss (int n		//!< order of aproximation
              ):_basis_order (n), 
                _alpha (_basis_order, _basis_order),
                _beta (_basis_order),_V(_basis_order,_basis_order),
                _S(_basis_order),_work(_basis_order),_x(_basis_order)
		{
      clear();
		};
    
    MNK_Gauss(const MNK_Gauss& a): _basis_order(a._basis_order), 
                                    _alpha (_basis_order, _basis_order),
                                    _beta (_basis_order),_V(_basis_order,_basis_order),
                                    _S(_basis_order),_work(_basis_order),_x(_basis_order)
    {
      clear();
    }
    
      
    const MNK_Gauss& operator=(const MNK_Gauss& a)
    {
      _basis_order=a._basis_order;
      _basis=a._basis;
      
      _alpha.resize(_basis_order,_basis_order);
      _beta.resize(_basis_order);
      _V.resize(_basis_order,_basis_order);
      _S.resize(_basis_order);
      _work.resize(_basis_order);
      _x.resize(_basis_order) ;
      clear();
    }
    
    void clear() 
    {
			_alpha = 0.0;
			_beta = 0.0;
      _V = 0.0;
      _S = 0.0;
      _x = 0.0;
      _work =0.0; // ?
    }
	
		//! fill the matrix, call as many times as there are points
		void accumulate(const T & pnt, double value)
		{
			for (int i = 0; i < dim (); i++) {
				for (int j = 0; j < dim (); j++) {
          _alpha.set(i,j, _alpha.get(i,j)+ _basis(i, pnt) * _basis(j, pnt));
				}
				//2: fill the right part
				_beta.set(i,_beta.get(i)+_basis(i, pnt) * value);
			}
		}
    
    //! fill the matrix, call as many times as there are points
    void accumulate_invw(const T & pnt, double value,double weight=1.0)
    {
      for (int i = 0; i < dim (); i++) {
        for (int j = 0; j < dim (); j++) {
          _alpha.set(i,j, _alpha.get(i,j)+ _basis(i, pnt) * _basis(j, pnt) * weight);
        }
        //2: fill the right part
        _beta.set(i,_beta.get(i)+_basis(i, pnt) * value * weight);
      }
    }
    
		//! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
		void accumulate(const valuesVec& basis, double value)
		{
			for (int i = 0; i < dim(); i++) {
				for (int j = 0; j < dim(); j++) {
          _alpha.set(i,j, _alpha.get(i,j)+basis[i] * basis[j]);
				}
				//2: fill the right part
				_beta.set(i,_beta.get(i)+basis[i] * value);
			}
		}
    
    //! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
    //! with inverted weights
    void accumulate_invw(const valuesVec& basis, double value,double w)
    {
      for (int i = 0; i < dim(); i++) {
        for (int j = 0; j < dim(); j++) {
          _alpha.set(i,j, _alpha.get(i,j)+basis[i] * basis[j]*w);
        }
        //2: fill the right part
        _beta.set(i,_beta.get(i)+basis[i] * value*w);
      }
    }

    
		//! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
		void accumulate(const valuesVecF& basis, double value)
		{
			for (int i = 0; i < dim(); i++) {
				for (int j = 0; j < dim(); j++) {
          _alpha.set(i,j, _alpha.get(i,j)+basis[i] * basis[j]);
				}
				//2: fill the right part
        _beta.set(i,_beta.get(i)+basis[i] * value);
			}
		}
	
		//! solve MNK aproximation problem - return vector of coeffecients
		int solve (valuesVec & coeff	//!< coeffecients
							)
		{
			//coeff.resize(dim ());	// make sure we have enough space!
			int i;
			for (i = 0; i < dim (); i++)
				coeff[i] = 0;
	
      //int r=gsl_linalg_SV_decomp(_alpha,_V,_S,_work);
      gsl_linalg_SV_decomp_jacobi(_alpha,_V,_S);
      gsl_linalg_SV_solve(_alpha,_V,_S,_beta,_x);
      //int r=gsl_linalg_QR_decomp(_alpha,_work);
      //gsl_linalg_QR_solve(_alpha,_work,_beta,_x);
      for (i = 0; i < dim (); i++)
				coeff[i] = _x.get(i);      
      return 0;
		}
    
    //! solve MNK aproximation problem - return vector of coeffecients
    int solve_unstable(valuesVec & coeff  //!< coeffecients
              )
    {
      //coeff.resize(dim ()); // make sure we have enough space!
      int i;
      for (i = 0; i < dim (); i++)
        coeff[i] = 0;
  
      gsl_linalg_SV_decomp(_alpha,_V,_S,_work);
      //gsl_linalg_SV_decomp_jacobi(_alpha,_V,_S);
      for (i = 0; i < dim (); i++)
        if(_S.get(i)<0.1) _S.set(i,0.0);
      
      gsl_linalg_SV_solve(_alpha,_V,_S,_beta,_x);
      //int r=gsl_linalg_QR_decomp(_alpha,_work);
      //gsl_linalg_QR_solve(_alpha,_work,_beta,_x);
      for (i = 0; i < dim (); i++)
        coeff[i] = _x.get(i);      
      return 0;
    }
    
  };

	template< class T > class Polinom
	{
	public:
		double operator() (int n, T x) const 
		{
			double res = 1.0;
			for (int i = 0; i < n; i++)
				res *= x;
			return res;
		}
	};
  
  //! modified version of polimolial, with no bias for order 1
	template< class T > class Polinom_Mod
	{
	public:
		double operator() (int n, T x) const
		{
      switch(n)
			{
        case 0: return x;
        case 1: return 1;
        default: //n>1
          double res = 1.0;
          for (int i = 0; i < n; i++)
            res *= x;
          return res;
      }
		}
	};
		
	typedef minc::MNK_Gauss<double, Polinom<double > >      MNK_Gauss_Polinomial;
  typedef minc::MNK_Gauss<double, Polinom_Mod<double > >  MNK_Gauss_Polinomial_Mod;
}; //minc

#endif //__GSL_GAUSS_H__
