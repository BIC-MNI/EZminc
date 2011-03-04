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
#ifndef __GSL_GLUE_H__
#define __GSL_GLUE_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <vector>

namespace minc
{
  //! gsl helper class
  class gsl_double_vector
  {
  protected:
    bool _own;
    void clean(void)
    {
      if(_vec && _own) gsl_vector_free(_vec);
        _vec=NULL;
    }
    
    void copy(const gsl_double_vector& a)
    {
      clean();
      if(a._vec)
      {
        resize(a._vec->size);
        gsl_vector_memcpy(_vec,a._vec);
      }
    }
 
    void copy(const std::vector<double>& a)
    {
      clean();
      if(!a.empty())
      {
        resize(a.size());
        for(int i=0;i<a.size();i++)
          gsl_vector_set(_vec,i,a[i]);
      }
    }
    
  public:
    gsl_vector *_vec;
    
    operator gsl_vector*()
    {
      return _vec;
    }
    
    operator const gsl_vector*() const
    {
      return _vec;
    }
    
    int size() const
    {
      if(!_vec) return 0;
      else return _vec->size;
    }
    gsl_double_vector():_vec(NULL)
    {}
    
    gsl_double_vector(gsl_vector* v):_vec(v),_own(false)
    {
    }
    
    void attach(gsl_vector* v)
    {
      clean();
      _vec=v;
      _own=false;
    }

    explicit gsl_double_vector(int n):_vec(NULL)
    {
      resize(n);
    }
    
    gsl_double_vector(int n,double v):_vec(NULL)
    {
      resize(n);
      set_all(v);
    }

    gsl_double_vector(const gsl_double_vector& a):_vec(NULL)
    {
      copy(a);
    }
    
    gsl_double_vector(const std::vector<double>& a):_vec(NULL)
    {
      copy(a);
    }
    
    const gsl_double_vector& operator=(const gsl_double_vector& a)
    {
      copy(a);
      return *this;
    }
    
    const gsl_double_vector& operator=(double a)
    {
      set_all(a);
      return *this;
    }
    
    const gsl_double_vector& operator=(const std::vector<double>& a)
    {
      copy(a);
      return *this;
    }
    
    const gsl_double_vector& operator*=(double a)
    {
      gsl_vector_scale(_vec,a);
      return *this;
    }
    
    
    ~gsl_double_vector()
    {
      clean();
    }
    void resize(int n)
    {
      clean();
      _vec=gsl_vector_alloc(n);
      _own=true;
    }
    
    double operator[](int n) const
    {
      return get(n);
    }
    
    double get(int n) const 
    {
      return gsl_vector_get(_vec,n);
    }
    
    double set(int n,double v) 
    {
      gsl_vector_set(_vec,n,v);
      return v;
    }
    
    void set_all(double v)
    {
      gsl_vector_set_all(_vec,v);
    }
    
    void copy_to(std::vector<double> &out) const
    {
      out.resize(size());
      for(int i=0;i<out.size();i++)
        out[i]=gsl_vector_get(_vec,i);
    }
    
    int read(FILE *f)
    {
      return gsl_vector_fread(f,_vec);
    }
    
    int write(FILE *f)
    {
      return gsl_vector_fwrite(f,_vec);
    }
  };
  
  
  //! GSL matrix
  class gsl_double_matrix
  {
  protected:
    bool _own;
  
    void clean(void)
    {
      if(_mtx && _own) gsl_matrix_free(_mtx);
        _mtx=NULL;
    }
    
    void copy(const gsl_double_matrix& a)
    {
      clean();
      if(a._mtx)
      {
        resize(a._mtx->size1,a._mtx->size2);
        gsl_matrix_memcpy(_mtx,a._mtx);
      }
    }
    
  public:
    gsl_matrix *_mtx;
  
    operator gsl_matrix *()
    {
      return _mtx;
    }
    
    operator const gsl_matrix *() const
    {
      return _mtx;
    }
  
    //!size of the ith dimension
    int size(int i) const
    {
      if(!_mtx) return 0;
      else return !i?_mtx->size1:_mtx->size2;
    }
    
    //! constructor
    gsl_double_matrix():_mtx(NULL)
    {
    }
    
    //! use another gsl matrix
    explicit gsl_double_matrix(gsl_matrix* m):_mtx(m),_own(false)
    {
    }
    
    //! attach another gsl matrix
    void attach(gsl_matrix* m)
    {
      clean();
      _mtx=m;
      _own=false;
    }
      
    //!  constructor, create matrix ixj
    explicit gsl_double_matrix(int i,int j):_mtx(NULL)
    {
      resize(i,j);
    }
    
    //!  constructor, create matrix ixj, set all elements to v
    gsl_double_matrix(int i,int j,double v):_mtx(NULL)
    {
      resize(i,j);
      set_all(v);
    }

    //! copy constructor
    gsl_double_matrix(const gsl_double_matrix& a):_mtx(NULL)
    {
      copy(a);
    }
    
    //! copy data from a std vector
    gsl_double_matrix(const std::vector<double>& a):_mtx(NULL)
    {
      copy(a);
    }
    
    //! copy values from another matrix
    const gsl_double_matrix& operator=(const gsl_double_matrix& a)
    {
      copy(a);
      return *this;
    }
    
    //! set all elements to a
    const gsl_double_matrix& operator=(double a)
    {
      set_all(a);
      return *this;
    }
    
    //! destructor
    ~gsl_double_matrix()
    {
      clean();
    }
    
    //! change matrix size, all data is lost
    void resize(int i,int j)
    {
      clean();
      _mtx=gsl_matrix_alloc(i,j);
      _own=true;
    }
    
    //! get element
    double get(int i,int j) const 
    {
      return gsl_matrix_get(_mtx,i,j);
    }
    
    //! set element
    double set(int i,int j,double v) 
    {
      gsl_matrix_set(_mtx,i,j,v);
      return v;
    }
    
    //! set all elements
    void set_all(double v)
    {
      gsl_matrix_set_all(_mtx,v);
    }
    
    //! return number of rows
    int rows(void) const
    {
      return _mtx->size1;
    }
    
    //! return number of columns
    int columns(void) const
    {
      return _mtx->size2;
    }
     
    //! read matrix from a file
    int read(FILE *f)
    {
      return gsl_matrix_fread(f,_mtx);
    }
    
    //! write matrix to a file
    int write(FILE *f)
    {
      return gsl_matrix_fwrite(f,_mtx);
    }
  };
  
  
  class _gsl_permutation
  {
    protected:
      gsl_permutation *_p;
      bool _own;
      
      void clean(void)
      {
        if(_p && _own) gsl_permutation_free(_p);
        _p=NULL;
      }
    
      void copy(const _gsl_permutation& a)
      {
        clean();
        if(a._p)
        {
          resize(a._p->size);
          gsl_permutation_memcpy(_p,a._p);
        }
      }
      
    public:
      
      ~_gsl_permutation()
      {
        clean();
      }
      
      _gsl_permutation():_p(NULL)
      {
      }
      
      explicit _gsl_permutation(gsl_permutation* p):_p(p),_own(false)
      {
      }
    
      void attach(gsl_permutation* p)
      {
        clean();
        _p=p;
        _own=false;
      }
      
      explicit _gsl_permutation(int i):_p(NULL)
      {
        resize(i);
      }

      _gsl_permutation(const _gsl_permutation& a):_p(NULL)
      {
        copy(a);
      }
      
      void resize(int i)
      {
        clean();
        _p=gsl_permutation_alloc(i);
        _own=true;
      }
      
      int size(void) const
      {
        if(_p) return _p->size;
        else return 0;
      }
      
      _gsl_permutation& operator=(const _gsl_permutation& a)
      {
        copy(a);
        return *this;
      }
      
      operator gsl_permutation* ()
      {
        return _p;
      }
      
      operator const gsl_permutation* () const
      {
        return _p;
      }
      
      

  };

  
}; //minc
#endif //__GSL_GLUE_H__
