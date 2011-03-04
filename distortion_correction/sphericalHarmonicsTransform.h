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
#ifndef __spherical_harmonics_transform_h__
#define __spherical_harmonics_transform_h__

#include <itkObject.h>
#include <itkPoint.h>
#include <itkVector.h>
#include <itkCovariantVector.h>
#include <vnl/vnl_vector_fixed.h>
#include <itkArray.h>
#include <itkArray2D.h>
#include <itkTransform.h>
#include <itkObjectFactory.h>
#include <memory>
#include "minc_wrappers.h"

#ifdef _DEBUG
#include <debug/vector>
#else
//#include <ext/pool_allocator.h>
#endif 

namespace minc
{

  typedef float basis_type;
  #ifdef _DEBUG
    typedef __gnu_debug::vector <basis_type> basis_vector;
  #else
    typedef std::vector<basis_type> basis_vector;
  //typedef std::vector<basis_type,__gnu_cxx::__pool_alloc<basis_type> > basis_vector;
  #endif
  
  class SphericalFunctions
  {
    public:
      double _scaling;
      SphericalFunctions():_scaling(100.0)
      {}
    
    double operator()(int n, tag_point p) const;
    static unsigned int parameters_no(int order);
    static double scale(int n,double);
    void generate_basis(basis_vector &basis, int order, tag_point p);
  };
  
  class CylindricalFunctions
  {
    public:
      double _scaling;
      CylindricalFunctions():_scaling(100.0)
      {}
    
    double operator()(int n, tag_point p) const;
    static unsigned int parameters_no(int order);
    static double scale(int n,double);
    void generate_basis(basis_vector &basis, int order, tag_point p);
  };
  
  template<class T> tag_point apply_transform(tag_point p,const basis_vector &basis, const T& parameters)
  {
    tag_point pnt;
    pnt.Fill(0.0);
    int _param_no=basis.size();
    for(int i=0;i<_param_no;i++)
    {
      double bs=basis[i];      
      pnt[0]+=bs*parameters[i];
      pnt[1]+=bs*parameters[i+_param_no];
      pnt[2]+=bs*parameters[i+_param_no*2];
    }
    return pnt;
  }
  
  template<class T> tag_point apply_transform3(tag_point p,const basis_vector &basis, const T& parameters)
  {
    tag_point pnt;
    pnt.Fill(0.0);
    int _param_no=basis.size();
    for(int i=0;i<_param_no;i++)
    {
      double bs=basis[i];      
      pnt[0]+=bs*parameters[0][i];
      pnt[1]+=bs*parameters[1][i];
      pnt[2]+=bs*parameters[2][i];
    }
    return pnt;
  }
  
  typedef SphericalFunctions basis_functions_x;
  typedef SphericalFunctions basis_functions_y;
  typedef SphericalFunctions basis_functions_z;


  /** \class SphericalHarmonicsTransform
   * \brief Implementation of an Spherical Harmonics Transform.
   *
   * \ingroup Transforms
   *
   */
  class SphericalHarmonicsTransform : public itk::Transform < double, 3, 3>
  {
  public:
    typedef double TScalarType;
    /** Standard class typedefs. */
    typedef SphericalHarmonicsTransform  Self;
    typedef itk::Transform< double, 3, 3 > Superclass;
    typedef itk::SmartPointer< Self >   Pointer;
    typedef itk::SmartPointer< const Self >  ConstPointer;
    
    
    /** New method for creating an object using a factory. */
    itkNewMacro(Self);
  
    /** Run-time type information (and related methods). */
    itkTypeMacro( SphericalHarmonicsTransform, itk::Transform );
  
    /** Dimension of the domain space. */
    itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
    
    /** Type of the input parameters. */
    
    typedef  double ScalarType;
  
    /** Type of the input parameters. */
    typedef Superclass::ParametersType ParametersType;
  
    /** Type of the Jacobian matrix. */
    typedef Superclass::JacobianType  JacobianType;
  
    /** Standard vector type for this class. */
    typedef itk::Vector<TScalarType,
                  itkGetStaticConstMacro(InputSpaceDimension)>  InputVectorType;
    typedef itk::Vector<TScalarType,
                  itkGetStaticConstMacro(OutputSpaceDimension)> OutputVectorType;
    
    /** Standard covariant vector type for this class */
    typedef itk::CovariantVector<TScalarType,
                            itkGetStaticConstMacro(InputSpaceDimension)>  InputCovariantVectorType;
                            
    typedef itk::CovariantVector<TScalarType,
                            itkGetStaticConstMacro(OutputSpaceDimension)> OutputCovariantVectorType;
    
    /** Standard coordinate point type for this class */
    typedef itk::Point<TScalarType,3 > InputPointType;
  
    typedef itk::Point<TScalarType,3 > OutputPointType;
    
    /**  Method to transform a point. */
    virtual OutputPointType TransformPoint(const InputPointType  &point ) const;
    OutputPointType TransformPointUnCached(const InputPointType  &point ) const;
  
    /**  Method to transform a vector. */
    
    virtual OutputVectorType TransformVector(const InputVectorType &vector) const 
    {
      itkExceptionMacro( << "Not Implemented" );
      return vector; 
    }
  
    /**  Method to transform a vnl_vector. */
    virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &vector) const
    {
      itkExceptionMacro( << "Not Implemented" );
      return vector; 
    }
  
    /**  Method to transform a CovariantVector. */
    virtual OutputCovariantVectorType TransformCovariantVector(
      const InputCovariantVectorType &vector) const
    {
      itkExceptionMacro( << "Not Implemented" );
      return vector; 
    }
  
    /** Set the transformation to an Identity
     */
    virtual void SetIdentity( void );
    void SetOrder(int order);
    
    virtual unsigned int GetNumberOfParameters(void) const
    {
      return _par_count*3;
    }
    
    /** Set the Transformation Parameters
     * and update the internal transformation. */
    virtual void  SetParameters(const ParametersType & param);
    virtual const ParametersType & GetParameters(void);
    
    void ImportParameters(const ParametersType & param,bool all=false);
    
    virtual const JacobianType & GetJacobian(const InputPointType  & point) const
    { 
      //this->m_Jacobian.Fill(0.0);
      minc::image3d::IndexType idx;
      _basis_cache->TransformPhysicalPointToIndex(point,idx);
      const basis_vector& bas=_basis_cache->GetPixel(idx);
      //if(bas.empty()) 
      //  {itkExceptionMacro( << "Basis not cached!" );}
      for(int i=0;i<_param_no;i++)
      {
        double bs=bas[i];
        this->m_Jacobian(0,i)=bs;
        this->m_Jacobian(1,i)=bs;
        this->m_Jacobian(2,i)=bs;
      }
      return this->m_Jacobian;
    }
    
  
    void calculate_basis(mask3d::Pointer sample);
    void calculate_basis(image3d::Pointer sample);
    void GetDeltas(ParametersType &delta);
    void GetScales(itk::Array< double > &scales);
    //void GetScales(basis_vector &scales);
    void SetParBaseCount(int b,int c);
    void SetCache(bool on=true)
    {
      _cache_on=on;
    }
    
    
  protected:
    SphericalHarmonicsTransform(): itk::Transform< double,3,3 >(3,3),
      _basis_cache(Basis_cache_vector::New()),_cache_on(false)
    { 
      _extent=100.0;
      _par_base=0;
      SetOrder(3);
      SetIdentity();
    }
    
    virtual ~SphericalHarmonicsTransform() {}
    
    ParametersType _parameters,_parameters2;
    
    SphericalFunctions basis;
    unsigned int _param_no;
    typedef itk::Image<basis_vector, 3 > Basis_cache_vector;
    typedef itk::ImageRegionIteratorWithIndex < Basis_cache_vector > basis_iterator; 
    mutable Basis_cache_vector::Pointer _basis_cache;
    mutable basis_vector _tmp;
    itk::Array< double > _scales;
    double _extent;
    bool _cache_on;
    int _par_base;
    int _par_count;
    
  private:
    SphericalHarmonicsTransform ( const Self & ); //purposely not implemented
    void operator= ( const Self & ); //purposely not implemented
  };

  inline void SphericalHarmonicsTransform::SetParameters(const SphericalHarmonicsTransform::ParametersType & param)
  {
    for(int i=0;i<_par_count;i++)
    {
      int j=i+_par_base;
       _parameters[j            ]=param[i];
       _parameters[j+_param_no  ]=param[i+_par_count];
       _parameters[j+_param_no*2]=param[i+_par_count*2];
    }
    
  }
  
  inline const SphericalHarmonicsTransform::ParametersType & SphericalHarmonicsTransform::GetParameters(void)
  {
    for(int i=0;i<_par_count;i++)
    {
      int j=i+_par_base;
       _parameters2[i             ]=_parameters[j];
       _parameters2[i+_par_count  ]=_parameters[j+_param_no];
       _parameters2[i+_par_count*2]=_parameters[j+_param_no*2];
    }
    return _parameters2;
  }
  
  
  inline SphericalHarmonicsTransform::OutputPointType SphericalHarmonicsTransform::TransformPoint(const SphericalHarmonicsTransform::InputPointType  &point ) const
  {
    OutputPointType pnt;
    pnt.Fill(0.0);
    minc::image3d::IndexType idx;
    _basis_cache->TransformPhysicalPointToIndex(point, idx);
    if(!_cache_on)
    {
      SphericalFunctions sph;
      sph.generate_basis(_tmp,_param_no,point);
      for(int i=0;i<_param_no;i++)
      {
        double bs=_tmp[i];
        pnt[0]+=bs*_parameters[i];
        pnt[1]+=bs*_parameters[i+_param_no];
        pnt[2]+=bs*_parameters[i+_param_no*2];
      }
    } else {
      basis_vector &bas= _basis_cache->GetPixel(idx);
      if(bas.empty())
        itkExceptionMacro( << "Trying to use unallocated cache element" );
      for(int i=0;i<_param_no;i++)
      {
        //TODO finish this
        double bs=bas[i];
        pnt[0]+=bs*_parameters[i];
        pnt[1]+=bs*_parameters[i+_param_no];
        pnt[2]+=bs*_parameters[i+_param_no*2];
      }
    }
    return pnt;
  }
  
  inline SphericalHarmonicsTransform::OutputPointType SphericalHarmonicsTransform::TransformPointUnCached(const SphericalHarmonicsTransform::InputPointType  &point ) const
  {
    OutputPointType pnt;
    pnt.Fill(0.0);
//    float tmp[200];
    SphericalFunctions sph;
    sph.generate_basis(_tmp,_param_no,point);
    for(int i=0;i<_param_no;i++)
    {
      double bs=_tmp[i];      
      pnt[0]+=bs*_parameters[i];
      pnt[1]+=bs*_parameters[i+_param_no];
      pnt[2]+=bs*_parameters[i+_param_no*2];
    }
    return pnt;
  }
  
  
  class CylindricalHarmonicsTransform : public itk::Transform < double, 3, 3>
  {
  public:
    typedef double TScalarType;
    /** Standard class typedefs. */
    typedef CylindricalHarmonicsTransform  Self;
    typedef itk::Transform< double, 3, 3 > Superclass;
    typedef itk::SmartPointer< Self >   Pointer;
    typedef itk::SmartPointer< const Self >  ConstPointer;
    
    
    /** New method for creating an object using a factory. */
    itkNewMacro(Self);
  
    /** Run-time type information (and related methods). */
    itkTypeMacro( SphericalHarmonicsTransform, itk::Transform );
  
    /** Dimension of the domain space. */
    itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
    
    /** Type of the input parameters. */
    
    typedef  double ScalarType;
  
    /** Type of the input parameters. */
    typedef Superclass::ParametersType ParametersType;
  
    /** Type of the Jacobian matrix. */
    typedef Superclass::JacobianType  JacobianType;
  
    /** Standard vector type for this class. */
    typedef itk::Vector<TScalarType,
                  itkGetStaticConstMacro(InputSpaceDimension)>  InputVectorType;
    typedef itk::Vector<TScalarType,
                  itkGetStaticConstMacro(OutputSpaceDimension)> OutputVectorType;
    
    /** Standard covariant vector type for this class */
    typedef itk::CovariantVector<TScalarType,
                            itkGetStaticConstMacro(InputSpaceDimension)>  InputCovariantVectorType;
                            
    typedef itk::CovariantVector<TScalarType,
                            itkGetStaticConstMacro(OutputSpaceDimension)> OutputCovariantVectorType;
    
    /** Standard coordinate point type for this class */
    typedef itk::Point<TScalarType,3 > InputPointType;
  
    typedef itk::Point<TScalarType,3 > OutputPointType;
    
    /**  Method to transform a point. */
    virtual OutputPointType TransformPoint(const InputPointType  &point ) const;
    OutputPointType TransformPointUnCached(const InputPointType  &point ) const;
  
    /**  Method to transform a vector. */
    
    virtual OutputVectorType TransformVector(const InputVectorType &vector) const 
    {
      itkExceptionMacro( << "Not Implemented" );
      return vector; 
    }
  
    /**  Method to transform a vnl_vector. */
    virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &vector) const
    {
      itkExceptionMacro( << "Not Implemented" );
      return vector; 
    }
  
    /**  Method to transform a CovariantVector. */
    virtual OutputCovariantVectorType TransformCovariantVector(
      const InputCovariantVectorType &vector) const
    {
      itkExceptionMacro( << "Not Implemented" );
      return vector; 
    }
  
    /** Set the transformation to an Identity
     */
    virtual void SetIdentity( void );
    void SetOrder(int order);
    
    virtual unsigned int GetNumberOfParameters(void) const
    {
      return _par_count*3;
    }
    
    /** Set the Transformation Parameters
     * and update the internal transformation. */
    virtual void  SetParameters(const ParametersType & param);
    virtual const ParametersType & GetParameters(void);
    
    void ImportParameters(const ParametersType & param,bool all=false);
    
    virtual const JacobianType & GetJacobian(const InputPointType  & point) const
    { 
      this->m_Jacobian.Fill(0.0);
      //TODO:finish this
      
      return this->m_Jacobian;
    }
    
  
    void calculate_basis(mask3d::Pointer sample);
    void calculate_basis(image3d::Pointer sample);
    void GetDeltas(ParametersType &delta);
    void GetScales(itk::Array< double > &scales);
    void SetParBaseCount(int b,int c);
    void SetCache(bool on=true)
    {
      _cache_on=on;
    }
    
    
  protected:
    CylindricalHarmonicsTransform(): itk::Transform< double,3,3 >(3,3),
      _basis_cache(Basis_cache_vector::New()),_cache_on(false)
    { 
      _extent=100.0;
      _par_base=0;
      SetOrder(3);
      SetIdentity();
    }
    
    virtual ~CylindricalHarmonicsTransform() {}
    
    ParametersType _parameters,_parameters2;
    
    CylindricalFunctions basis;
    unsigned int _param_no;
    typedef itk::Image<basis_vector, 3 > Basis_cache_vector;
    typedef itk::ImageRegionIteratorWithIndex < Basis_cache_vector > basis_iterator; 
    mutable Basis_cache_vector::Pointer _basis_cache;
    mutable basis_vector _tmp;
    itk::Array< double > _scales;
    double _extent;
    bool _cache_on;
    int _par_base;
    int _par_count;
    
  private:
    CylindricalHarmonicsTransform ( const Self & ); //purposely not implemented
    void operator= ( const Self & ); //purposely not implemented
  };

  inline void CylindricalHarmonicsTransform::SetParameters(const CylindricalHarmonicsTransform::ParametersType & param)
  {
    for(int i=0;i<_par_count;i++)
    {
      int j=i+_par_base;
       _parameters[j            ]=param[i];
       _parameters[j+_param_no  ]=param[i+_par_count];
    }
    
  }
  
  inline const CylindricalHarmonicsTransform::ParametersType & CylindricalHarmonicsTransform::GetParameters(void)
  {
    for(int i=0;i<_par_count;i++)
    {
      int j=i+_par_base;
       _parameters2[i             ]=_parameters[j];
       _parameters2[i+_par_count  ]=_parameters[j+_param_no];
    }
    return _parameters2;
  }
  
  
  inline CylindricalHarmonicsTransform::OutputPointType CylindricalHarmonicsTransform::TransformPoint(const CylindricalHarmonicsTransform::InputPointType  &point ) const
  {
    OutputPointType pnt;
    pnt.Fill(0.0);
    minc::image3d::IndexType idx;
    _basis_cache->TransformPhysicalPointToIndex(point, idx);
    double rp=sqrt(point[0]*point[0]+point[1]*point[1]);
    if(!_cache_on)
    {
      CylindricalFunctions sph;
      sph.generate_basis(_tmp,_param_no,point);
      double r=0;
      for(int i=0;i<_param_no;i++)
      {
        double bs=_tmp[i];
        r+=bs*_parameters[i];
        pnt[2]+=bs*_parameters[i+_param_no];
      }
      if(rp>1e-6)
      {
        pnt[0]=point[0]*r/rp;
        pnt[1]=point[1]*r/rp;
      }else{
        pnt[0]=point[0];
        pnt[1]=point[1];
      }
    } else {
      basis_vector &bas= _basis_cache->GetPixel(idx);
      if(bas.empty())
        itkExceptionMacro( << "Trying to use unallocated cache element" );
      double r=0;
      for(int i=0;i<_param_no;i++)
      {
        double bs=bas[i];
        r+=bs*_parameters[i];
        pnt[2]+=bs*_parameters[i+_param_no];
      }
      if(rp>1e-6)
      {
        pnt[0]=point[0]*r/rp;
        pnt[1]=point[1]*r/rp;
      }else{
        pnt[0]=point[0];
        pnt[1]=point[1];
      }
    }
    return pnt;
  }
  
  inline CylindricalHarmonicsTransform::OutputPointType 
    CylindricalHarmonicsTransform::TransformPointUnCached(const CylindricalHarmonicsTransform::InputPointType  &point ) const
  {
    OutputPointType pnt;
    pnt.Fill(0.0);
//    float tmp[200];
    CylindricalFunctions sph;
    sph.generate_basis(_tmp,_param_no,point);
    double r=0;
    double rp=sqrt(point[0]*point[0]+point[1]*point[1]);
    for(int i=0;i<_param_no;i++)
    {
      double bs=_tmp[i];
      r+=bs*_parameters[i];
      pnt[2]+=bs*_parameters[i+_param_no];
    }
    if(rp>1e-6)
    {
      pnt[0]=point[0]*r/rp;
      pnt[1]=point[1]*r/rp;
    }else{
      pnt[0]=point[0];
      pnt[1]=point[1];
    }
    return pnt;
  }
  
}; // end namespace minc


#endif //__spherical_harmonics_transform_h__



