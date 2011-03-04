/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2006 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#ifndef __DATA_PROC_H__
#define __DATA_PROC_H_

#include <algorithm>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>

/*#include <vnl/algo/vnl_powell.h>
#include <vnl/algo/vnl_amoeba.h>
#include <vnl/algo/vnl_levenberg_marquardt.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_least_squares_function.h>
#include <vnl/vnl_least_squares_cost_function.h>*/

#include "minc_wrappers.h"

namespace minc
{
	template < class T, class B > class MNK_Gauss
	{
	protected:
		//basis_function basis;
		int _basis_order;
		vnl_matrix < double > _alpha;
		vnl_vector < double > _beta;
		B _basis;
		typedef std::vector < double > valuesVec;
    typedef std::vector < float > valuesVecF;
	public:
		double fit (const valuesVec & coeff, T pnt)
		{
			double val = 0.0;
			for (int i = 0; i < dim (); i++)
				val += _basis (i, pnt) * coeff[i];
			return val;
		}
		double fit (const valuesVec& basis,const valuesVec & coeff, T pnt)
		{
			double val = 0.0;
			for (int i = 0; i < dim (); i++)
				val += basis[i] * coeff[i];
			return val;
		}
	
		double fit (const valuesVec & coeff, T pnt, T min, T max)
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
							):_basis_order (n), _basis (basis), _alpha (n, n), _beta (n)
		{
			_alpha = 0.0;
			_beta = 0.0;
		};
	
		//! basic construcor
		MNK_Gauss (int n		//!< order of aproximation
							): _basis_order (n), _alpha (n, n), _beta (n)
		{
			_alpha = 0.0;
			_beta = 0.0;
		};
    
    MNK_Gauss(const MNK_Gauss& a):_basis_order(a._basis_order), _alpha(a._basis_order, a._basis_order), _beta(a._basis_order)
    {
			_alpha = 0.0;
			_beta = 0.0;
    }
      
    const MNK_Gauss& operator=(const MNK_Gauss& a)
    {
      _basis_order=a._basis_order;
      _alpha.set_size(_basis_order,_basis_order);
      _beta.set_size(_basis_order);
			_alpha = 0.0;
			_beta = 0.0;
    }
    
		void clear (void)
		{
			_alpha = 0.0;
			_beta = 0.0;
		}
	
		//! fill the matrix, call as many times as there are points
		void accumulate(const T & pnt, double value)
		{
			for (int i = 0; i < dim (); i++) {
				for (int j = 0; j < dim (); j++) {
					_alpha[i][j] += _basis(i, pnt) * _basis(j, pnt);
				}
				//2: fill the right part
				_beta[i] += _basis(i, pnt) * value;
			}
		}
    
		//! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
		void accumulate(const valuesVec& basis, double value)
		{
			for (int i = 0; i < dim(); i++) {
				for (int j = 0; j < dim(); j++) {
					_alpha[i][j] += basis[i] * basis[j];
				}
				//2: fill the right part
				_beta[i] += basis[i] * value;
			}
		}

    
		//! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
		void accumulate(const valuesVecF& basis, double value)
		{
			for (int i = 0; i < dim(); i++) {
				for (int j = 0; j < dim(); j++) {
					_alpha[i][j] += basis[i] * basis[j];
				}
				//2: fill the right part
				_beta[i] += basis[i] * value;
			}
		}
	
		//! solve MNK aproximation problem - return vector of coeffecients
		int solve (valuesVec & coeff	//!< coeffecients
							)
		{
			coeff.resize (dim ());	// make sure we have enough space!
			int i;
			for (i = 0; i < dim (); i++)
				coeff[i] = 0;
	
			//vnl_svd< double > svd(_alpha);
      vnl_qr< double > eq(_alpha);
			vnl_vector< double > solution(dim());
			solution = eq.solve(_beta);
			for (i = 0; i < dim (); i++)
				coeff[i] = solution[i];
			return coeff.size();//svd.rank ();
		};
    
    
		//! solve MNK aproximation problem - return vector of coeffecients
		double solve_svd(valuesVec & coeff	//!< coeffecients
							)
		{
			coeff.resize (dim ());	// make sure we have enough space!
			int i;
			for (i = 0; i < dim (); i++)
				coeff[i] = 0;
	
			vnl_svd< double > svd(_alpha);
      //vnl_qr< double > eq(_alpha);
			vnl_vector< double > solution(dim());
			solution = svd.solve(_beta);
			for (i = 0; i < dim (); i++)
				coeff[i] = solution[i];
			return svd.sigma_max()/svd.sigma_min();
		};
    
	};

	template < class T> class MNK_Gauss_opt
	{
	protected:
		//basis_function basis;
		int _basis_order;
		vnl_matrix < double > _alpha;
		vnl_vector < double > _beta;
		typedef std::vector < double > valuesVec;
    typedef std::vector < float > valuesVecF;
	public:
    
		double fit (const valuesVec& basis,const valuesVec & coeff, T pnt)
		{
			double val = 0.0;
			for (int i = 0; i < dim (); i++)
				val += basis[i] * coeff[i];
			return val;
		}
	
		int dim (void) const
		{
			return _basis_order;
		}
	
		//! basic construcor
		MNK_Gauss_opt (int n		//!< order of aproximation
							): _basis_order (n), _alpha (n, n), _beta (n)
		{
			_alpha = 0.0;
			_beta = 0.0;
		};
    
    MNK_Gauss_opt(const MNK_Gauss_opt& a):_basis_order(a._basis_order), _alpha(a._basis_order, a._basis_order), _beta(a._basis_order)
    {
			_alpha = 0.0;
			_beta = 0.0;
    }
      
    const MNK_Gauss_opt& operator=(const MNK_Gauss_opt& a)
    {
      _basis_order=a._basis_order;
      _alpha.set_size(_basis_order,_basis_order);
      _beta.set_size(_basis_order);
			_alpha = 0.0;
			_beta = 0.0;
    }
    
		void clear (void)
		{
			_alpha = 0.0;
			_beta = 0.0;
		}
	    
		//! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
		void accumulate(const double* basis, double value)
		{
			for (int i = 0; i < dim(); i++) {
				for (int j = 0; j < dim(); j++) {
					_alpha[i][j] += basis[i] * basis[j];
				}
				//2: fill the right part
				_beta[i] += basis[i] * value;
			}
		}

		//! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
		void accumulate(const valuesVec& basis, double value)
		{
			for (int i = 0; i < dim(); i++) {
				for (int j = 0; j < dim(); j++) {
					_alpha[i][j] += basis[i] * basis[j];
				}
				//2: fill the right part
				_beta[i] += basis[i] * value;
			}
		}
    
		//! fill the matrix, call as many times as there are points
    //! basis functions are precalculated
		void accumulate(const valuesVecF& basis, double value)
		{
			for (int i = 0; i < dim(); i++) {
				for (int j = 0; j < dim(); j++) {
					_alpha[i][j] += basis[i] * basis[j];
				}
				//2: fill the right part
				_beta[i] += basis[i] * value;
			}
		}
	
		//! solve MNK aproximation problem - return vector of coeffecients
		int solve (valuesVec & coeff	//!< coeffecients
							)
		{
			coeff.resize (dim ());	// make sure we have enough space!
			int i;
			for (i = 0; i < dim (); i++)
				coeff[i] = 0;
	
			//vnl_svd< double > svd(_alpha);
      vnl_qr< double > eq(_alpha);
			vnl_vector< double > solution(dim());
			solution = eq.solve(_beta);
			for (i = 0; i < dim (); i++)
				coeff[i] = solution[i];
			return coeff.size();//svd.rank ();
		};
    
		//! solve MNK aproximation problem - return vector of coeffecients
		int solve (double * coeff	//!< coeffecients
							)
		{
			//coeff.resize (dim ());	// make sure we have enough space!
			int i;
			for (i = 0; i < dim (); i++)
				coeff[i] = 0;
	
			//vnl_svd< double > svd(_alpha);
      vnl_qr< double > eq(_alpha);
			vnl_vector< double > solution(dim());
			solution = eq.solve(_beta);
			for (i = 0; i < dim (); i++)
				coeff[i] = solution[i];
			return dim();//svd.rank ();
		};
    
		//! solve MNK aproximation problem - return vector of coeffecients
		double solve_svd(double * coeff	//!< coeffecients
							)
		{
			//coeff.resize (dim ());	// make sure we have enough space!
			int i;
			for (i = 0; i < dim (); i++)
				coeff[i] = 0;
	
			vnl_svd< double > svd(_alpha);
      //vnl_qr< double > eq(_alpha);
			vnl_vector< double > solution(dim());
			solution = svd.solve(_beta);
			for (i = 0; i < dim (); i++)
				coeff[i] = solution[i];
			return svd.sigma_max()/svd.sigma_min();
		};
    
		//! solve MNK aproximation problem - return vector of coeffecients
		double solve_svd(valuesVec & coeff	//!< coeffecients
							)
		{
			coeff.resize (dim ());	// make sure we have enough space!
			int i;
			for (i = 0; i < dim (); i++)
				coeff[i] = 0;
	
			vnl_svd< double > svd(_alpha);
      //vnl_qr< double > eq(_alpha);
			vnl_vector< double > solution(dim());
			solution = svd.solve(_beta);
			for (i = 0; i < dim (); i++)
				coeff[i] = solution[i];
			return svd.sigma_max()/svd.sigma_min();
		};
	};

	template< class T > class Polinom
	{
	public:
		double operator() (int n, T x)
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
		double operator() (int n, T x)
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

		
	typedef minc::MNK_Gauss<double, Polinom<double > >  MNK_Gauss_Polinomial;
  typedef minc::MNK_Gauss<double, Polinom_Mod<double > >  MNK_Gauss_Polinomial_Mod;

  class General2ndOrderdPolinomial
  {
  public:
    double operator() (int n, tag_point p)
    {
      switch(n)
      {
        default:
        case 0: return 1;
        case 1: return p[0];
        case 2: return p[1];
        case 3: return p[2];
        case 4: return p[0]*p[0];
        case 5: return p[1]*p[1];
        case 6: return p[2]*p[2];
        case 7: return p[0]*p[1];
        case 8: return p[1]*p[2];
        case 9: return p[2]*p[0];
      }
    }
  };
  
  typedef minc::MNK_Gauss < tag_point, General2ndOrderdPolinomial >  MNK_Gauss_Polinomial3D;
	
};//minc

#endif //__DATA_PROC_H_
