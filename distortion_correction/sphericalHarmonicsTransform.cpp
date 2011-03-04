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
#include "sphericalHarmonicsTransform.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_result.h>
#include <iostream>

using namespace std;

namespace minc
{
  
	inline double SQR2(double x)
	{
		return x*x;
	}

	inline double Bz_generate(int n,int m,int part,double x,double y,double z)
	{
		double rxy = sqrt(SQR2(x) + SQR2(y));
		double r =   sqrt(SQR2(x) + SQR2(y) + SQR2(z));
		const double epsilon=1e-5;
		double costheta=r>epsilon ? z/r : 1.0 ; //z?
		
		double cosphi1,sinphi1,cosphi,sinphi,cosphi_,sinphi_;
		double rr=1.0;
		int i;
		if(!m) {
			cosphi1=cosphi=1.0;
			sinphi1=sinphi=0.0;
		} else {
			if(rxy>epsilon) {
				cosphi1=cosphi=x/rxy; 
				sinphi1=sinphi=y/rxy;
			} else {
				cosphi1=cosphi=1.0;
				sinphi1=sinphi=0.0;
			}
		}
		for(i = 2; i <= m; i++) {
			cosphi_ = (cosphi1 * cosphi) - (sinphi1 * sinphi);
			sinphi_ = (cosphi1 * sinphi) + (sinphi1 * cosphi);
			cosphi=cosphi_;
			sinphi=sinphi_;
		}
		for(i=1;i<=n;i++) {
			rr*=r;
		}
		gsl_sf_result tmp;
    //std::cerr << "gsl_sf_legendre_sphPlm_e("<<n << "," << m << "," << costheta << ")" << std::endl;
		if(gsl_sf_legendre_sphPlm_e(n, m, costheta, &tmp) != GSL_SUCCESS)
    {
      REPORT_ERROR("Error in gsl_sf_legendre_sphPlm_e");
    }
		if(part)
			return rr*tmp.val*sinphi;
		else
			return rr*tmp.val*cosphi;
	}
  
	
	double SphericalFunctions::operator() (int n, tag_point p) const
	{
    REPORT_ERROR("Intentionally not implemented!");
	}
  
	void SphericalFunctions::generate_basis(basis_vector &basis, int order, tag_point p)
  //this could be heavely optimized!
  {
    basis.resize(order);
    
		double rxy = sqrt(SQR2(p[0]) + SQR2(p[1]));
		double r =   sqrt(SQR2(p[0]) + SQR2(p[1]) + SQR2(p[2]));
		const double epsilon= 1e-6;
		double costheta= r > epsilon ? p[2]/r : 1.0 ; //z?
		
		double cosphi1,sinphi1,cosphi,sinphi,cosphi_,sinphi_;
		double rr=1.0;
    if(rxy>epsilon) {
      cosphi1=p[0]/rxy; 
      sinphi1=p[1]/rxy;
    } else {
      cosphi1=1.0;
      sinphi1=0.0;
    }

    //double sc=100.0;
    r/=_scaling;
    int n=2,j=0,i=0;
    basis[0]=p[2];
    basis[1]=p[0];
    basis[2]=p[1];
    cosphi=1.0;
    sinphi=0.0;
    rr=r*r*_scaling;
    for(int c=3;c<order;)
    {
      gsl_sf_result tmp;
      if(gsl_sf_legendre_sphPlm_e(n, j, costheta, &tmp) != GSL_SUCCESS)
        REPORT_ERROR("Error in gsl_sf_legendre_sphPlm_e");
      basis[c]=  rr*tmp.val*cosphi;
      c++;
      if(j)
      {
       if(c>=order)
         REPORT_ERROR("Out of bounds!");
       basis[c]=rr*tmp.val*sinphi;
       c++;
      }
      j++;
      if(j==(n+1)) 
      { 
        j=0; n++;
        rr*=r;
        cosphi=1.0;
        sinphi=0.0;
      } else {
        cosphi_ = (cosphi1 * cosphi) - (sinphi1 * sinphi);
        sinphi_ = (cosphi1 * sinphi) + (sinphi1 * cosphi);
        cosphi= cosphi_;
        sinphi= sinphi_;
      }
    }
  }
  
  double SphericalFunctions::scale(int n, double v)
	{
    if(n<3) return v;
    return 
      pow(v,floor(-1.0+sqrt(1.0+n)));
	}
  
	unsigned int SphericalFunctions::parameters_no(int order)
	{ 
    return order*order+2*order;
	}

  
	void CylindricalFunctions::generate_basis(basis_vector &basis, int order, tag_point p)
  //this could be heavely optimized!
  {
    basis.resize(order);
    
		double rxy = sqrt(SQR2(p[0]) + SQR2(p[1]));
		double r =   sqrt(SQR2(p[0]) + SQR2(p[1]) + SQR2(p[2]));
		const double epsilon= 1e-6;
		double costheta= r > epsilon ? p[2]/r : 1.0 ; //z?
		
		double cosphi1,sinphi1,cosphi,sinphi,cosphi_,sinphi_;
		double rr=1.0;
    //double sc=100.0;
    r/=_scaling;
    int n=2,j=0,i=0;
    
    basis[0]=rxy;
    basis[1]=p[2];
    rr=r*r*_scaling;
    for(int c=2;c<order;c++)
    {
      gsl_sf_result tmp;
      if(gsl_sf_legendre_sphPlm_e(n, j, costheta, &tmp) != GSL_SUCCESS)
        REPORT_ERROR("Error in gsl_sf_legendre_sphPlm_e");
      basis[c]=  rr*tmp.val;
      j++;
      if(j==(n+1)) 
      { 
        j=0; n++;
        rr*=r;
      }
    }
  }
  
  double CylindricalFunctions::scale(int n, double v)
	{
    if(n<2) return v;
    return 
      pow(v,floor(-0.75+sqrt(1.5*1.5+8*n)/2));
	}
  
	unsigned int CylindricalFunctions::parameters_no(int order)
	{ 
    return order*(order+1)/2+order;
	}
  
  void SphericalHarmonicsTransform::ImportParameters(const ParametersType & param, bool all)
  {
    int pno=param.Size()/3;

    if(pno>_param_no) // extend???
    {
      if(all)
      {
       _param_no=pno;
       _par_count=_param_no;
       _parameters.SetSize(_param_no*3);
       _parameters2.SetSize(_par_count*3);
      } else {
        pno=_param_no;
      }
    }
    _parameters.Fill(0.0);
		for(int i=0;i<pno;i++)
		{
			_parameters[i            ]=param[i];
			_parameters[i+_param_no  ]=param[i+pno];
			_parameters[i+_param_no*2]=param[i+pno*2];;
		}
  }
  
  void SphericalHarmonicsTransform::SetIdentity( void )
  {
    _parameters.Fill(0.0);
		_parameters[1          ]=1.0; //X
		_parameters[_param_no+2]=1.0; //Y
		_parameters[_param_no*2]=1.0; //Z
  }

  void SphericalHarmonicsTransform::SetOrder(int order)
  {
    _param_no=basis.parameters_no( order );
    _par_base=0;
    _par_count=_param_no;
    _parameters.SetSize( _param_no*3 );
    SetIdentity();
    this->m_Jacobian = JacobianType(3,_param_no);
    _tmp.resize( _param_no );
    //_scales.resize(_param_no);
    GetScales( _scales );
    _parameters.SetSize(_par_count*3);
    _parameters2.SetSize(_par_count*3);
    
  }

	void SphericalHarmonicsTransform::calculate_basis(mask3d::Pointer sample)
	{
    allocate_same( _basis_cache, sample );
		basis_iterator  itb(_basis_cache , sample->GetRequestedRegion() );
		mask3d_iterator it ( sample      , sample->GetRequestedRegion() );
    _cache_on=true;
    
		for(it.GoToBegin();!it.IsAtEnd();++it,++itb)
		{
			//if( !it.Value() ) 
      //  continue;
      //if( !itb.Value().empty() ) continue; // assume that we have already calculated coefficients here
      image3d::IndexType idx=it.GetIndex();
      tag_point p;
      sample->TransformIndexToPhysicalPoint(idx, p);
      SphericalFunctions sph;
      sph.generate_basis(_tmp, _param_no, p);
      itb.Set(_tmp);
		}
	}

	void SphericalHarmonicsTransform::calculate_basis(image3d::Pointer sample)
	{
    allocate_same(_basis_cache,sample);
    _cache_on=true;
		basis_iterator   itb(_basis_cache,sample->GetRequestedRegion());
		image3d_iterator it(sample,sample->GetRequestedRegion());
		for(it.GoToBegin();!it.IsAtEnd();++it,++itb)
		{
      image3d::IndexType idx=it.GetIndex();
      if( !itb.Value().empty() ) continue; // assume that we have already      
      tag_point p;
      sample->TransformIndexToPhysicalPoint(idx,p);
      SphericalFunctions sph;
      sph.generate_basis(_tmp,_param_no,p);
      itb.Set(_tmp);
#ifdef _DEBUG
#endif //_DEBUG
		}
	}

	void SphericalHarmonicsTransform::GetDeltas(ParametersType &delta)
	{
		delta.SetSize(GetNumberOfParameters());
		for(int i=0;i<_par_count;i++)
		{
			delta[i             ]=0.01; 
			delta[i+_par_count  ]=0.01;
			delta[i+_par_count*2]=0.01; 
		}
	}
  
	void SphericalHarmonicsTransform::GetScales(itk::Array< double > &scales)
	{
		scales.SetSize(GetNumberOfParameters());
		for(int i=0;i<_par_count;i++)
		{
			scales[i            ]= _extent/SphericalFunctions::scale(i+_par_base, _extent);
			scales[i+_par_count  ]= _extent/SphericalFunctions::scale(i+_par_base, _extent);
			scales[i+_par_count*2]= _extent/SphericalFunctions::scale(i+_par_base, _extent);
		}
	}
  
	void SphericalHarmonicsTransform::SetParBaseCount(int b, int c)
  {
    _par_base=b;
    if(c>0) _par_count=c;
    if(c<=0 || _par_count>(_param_no-_par_base)) _par_count=_param_no-_par_base;
    _parameters2.SetSize(_par_count*3);
  }

// ----
  void CylindricalHarmonicsTransform::ImportParameters(const ParametersType & param, bool all)
  {
    int pno=param.Size()/2;

    if(pno>_param_no) // extend???
    {
      if(all)
      {
       _param_no=pno;
       _par_count=_param_no;
       _parameters.SetSize(_param_no*2);
      } else {
        pno=_param_no;
      }
    }
    _parameters.Fill(0.0);
		for(int i=0;i<pno;i++)
		{
			_parameters[i            ]=param[i];
			_parameters[i+_param_no  ]=param[i+pno];
		}
  }
  
  void CylindricalHarmonicsTransform::SetIdentity( void )
  {
    _parameters.Fill(0.0);
		_parameters[0          ]=1.0; //R
		_parameters[_param_no+1]=1.0; //Z
  }

  void CylindricalHarmonicsTransform::SetOrder(int order)
  {
    _param_no=basis.parameters_no( order );
    _par_base=0;
    _par_count=_param_no;
    _parameters.SetSize( _param_no*2 );
    SetIdentity();
    this->m_Jacobian = JacobianType(3,_param_no);
    _tmp.resize( _param_no );
    //_scales.resize(_param_no);
    GetScales( _scales );
    _parameters.SetSize(_par_count*2);
    _parameters2.SetSize(_par_count*2);
    
  }

	void CylindricalHarmonicsTransform::calculate_basis(mask3d::Pointer sample)
	{
    allocate_same( _basis_cache, sample );
		basis_iterator  itb(_basis_cache , sample->GetRequestedRegion() );
		mask3d_iterator it ( sample      , sample->GetRequestedRegion() );
    _cache_on=true;
    
		for(it.GoToBegin();!it.IsAtEnd();++it,++itb)
		{
			//if( !it.Value() ) 
      //  continue;
      //if( !itb.Value().empty() ) continue; // assume that we have already calculated coefficients here
      image3d::IndexType idx=it.GetIndex();
      tag_point p;
      sample->TransformIndexToPhysicalPoint(idx, p);
      SphericalFunctions sph;
      sph.generate_basis(_tmp, _param_no, p);
      itb.Set(_tmp);
		}
	}

	void CylindricalHarmonicsTransform::calculate_basis(image3d::Pointer sample)
	{
    allocate_same(_basis_cache,sample);
    _cache_on=true;
		basis_iterator   itb(_basis_cache,sample->GetRequestedRegion());
		image3d_iterator it(sample,sample->GetRequestedRegion());
		for(it.GoToBegin();!it.IsAtEnd();++it,++itb)
		{
      image3d::IndexType idx=it.GetIndex();
      if( !itb.Value().empty() ) continue; // assume that we have already      
      tag_point p;
      sample->TransformIndexToPhysicalPoint(idx,p);
      CylindricalFunctions sph;
      sph.generate_basis(_tmp,_param_no,p);
      itb.Set(_tmp);
#ifdef _DEBUG
#endif //_DEBUG
		}
	}

	void CylindricalHarmonicsTransform::GetDeltas(ParametersType &delta)
	{
		delta.SetSize(GetNumberOfParameters());
		for(int i=0;i<_par_count;i++)
		{
			delta[i             ]=0.01; 
			delta[i+_par_count  ]=0.01;
		}
	}
  
	void CylindricalHarmonicsTransform::GetScales(itk::Array< double > &scales)
	{
		scales.SetSize(GetNumberOfParameters());
		for(int i=0;i<_par_count;i++)
		{
			scales[i           ]= _extent/CylindricalFunctions::scale(i+_par_base, _extent);
			scales[i+_par_count]= _extent/CylindricalFunctions::scale(i+_par_base, _extent);
		}
	}
  
	void CylindricalHarmonicsTransform::SetParBaseCount(int b, int c)
  {
    _par_base=b;
    if(c>0) _par_count=c;
    if(c<=0 || _par_count>(_param_no-_par_base)) _par_count=_param_no-_par_base;
    _parameters2.SetSize(_par_count*3);
  }
};
