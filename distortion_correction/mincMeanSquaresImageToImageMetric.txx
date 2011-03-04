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

#ifndef _mincMeanSquaresImageToImageMetric_txx
#define _mincMeanSquaresImageToImageMetric_txx

#include "mincMeanSquaresImageToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
//#include "mt19937ar.h"

namespace minc
{

	/*
	 * Constructor
	 */
	template <class TFixedImage, class TMovingImage,class TFixedTemplate> 
	MeanSquaresImageToImageMetric<TFixedImage,TMovingImage,TFixedTemplate>
	::MeanSquaresImageToImageMetric()
	{
    _keep=1.0;
		itkDebugMacro("Constructor");
		_template=0;
    _rng=gsl_rng_alloc(gsl_rng_default);
	}
  
  
	/*
  * Constructor
  */
  template <class TFixedImage, class TMovingImage,class TFixedTemplate> 
      MeanSquaresImageToImageMetric<TFixedImage,TMovingImage,TFixedTemplate>
  ::~MeanSquaresImageToImageMetric()
  {
    if(_rng) gsl_rng_free(_rng);
  }
  


	template <class TFixedImage, class TMovingImage, class TFixedTemplate> 
	typename MeanSquaresImageToImageMetric< TFixedImage, TMovingImage, TFixedTemplate >::MeasureType MeanSquaresImageToImageMetric<TFixedImage,TMovingImage,TFixedTemplate> ::GetMean( const TransformParametersType & parameters, const  FixedImageType* fixed, const  InterpolatorType* interpolator, const FixedTemplateType* mask) const
	{
		itkDebugMacro("GetMean( " << parameters<<" , "<< fixed <<" , "<< interpolator << " , "<< mask << " ) ");
		
		typedef  itk::ImageRegionConstIteratorWithIndex<TFixedTemplate> MaskIteratorType;
		
		MaskIteratorType ti( mask, this->GetFixedImageRegion() );
		
    typename FixedImageType::IndexType index;
	
		MeasureType measure = itk::NumericTraits< MeasureType >::Zero;
	
		int cnt = 0;
		
		this->SetTransformParameters( parameters );
	
		for(ti.GoToBegin();!ti.IsAtEnd();++ti)
		{
      if(ti.Get() && (_keep>=1.0 || gsl_rng_uniform(_rng)<_keep))
			{
				index = ti.GetIndex();
				
				typename Superclass::InputPointType inputPoint;
				fixed->TransformIndexToPhysicalPoint( index, inputPoint );
		
				typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );
		
				if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
				{
					const RealType movingValue  = interpolator->Evaluate( transformedPoint );
					const RealType fixedValue =   fixed->GetPixel(index);
					cnt++;
					const RealType diff = movingValue - fixedValue; 
					measure += diff * diff; 
				}
			}
		}
		
		if( !cnt)
		{
			//itkExceptionMacro(<<"All the points mapped to outside of the moving image");
      //return penalty value!
      return 10;
		}
    
		else
		{
			measure /= cnt;
		}
		return measure;
		
	}

	
	/*
	 * Get the match Measure
	 */
	template <class TFixedImage, class TMovingImage, class TFixedTemplate> 
	typename MeanSquaresImageToImageMetric< TFixedImage, TMovingImage, TFixedTemplate >::MeasureType MeanSquaresImageToImageMetric<TFixedImage,TMovingImage,TFixedTemplate> ::GetValue( const TransformParametersType & parameters ) const
	{
		itkDebugMacro("GetValue( " << parameters << " ) ");
    MeasureType m= GetMean(parameters, this->m_FixedImage, this->m_Interpolator, _template);
    
    if(_template2)
    {
      m+=GetMean(parameters, _fixed2, _interpolator2, _template2);
    }
    return m;
	}
	
	/*
	 * Get the Derivative Measure
	 */
	template < class TFixedImage, class TMovingImage, class TFixedTemplate> 
	void
	MeanSquaresImageToImageMetric<TFixedImage,TMovingImage,TFixedTemplate>
	::GetDerivative( const TransformParametersType & parameters,
									 DerivativeType & derivative  ) const
	{
		itkExceptionMacro( << "Not Implemented" );
	  //std::cout<<"+"<<std::flush;
		itkDebugMacro("GetDerivative( " << parameters << " ) ");
		
		if( !this->GetGradientImage() )
		{
			itkExceptionMacro(<<"The gradient image is null, maybe you forgot to call Initialize()");
		}
	
		FixedImageConstPointer fixedImage = this->m_FixedImage;
	
		if( !fixedImage || !_template ) 
    {
			itkExceptionMacro( << "Fixed image has not been assigned" );
		}
	
		const unsigned int ImageDimension = FixedImageType::ImageDimension;
	
	
		typedef itk::ImageRegionConstIteratorWithIndex<FixedTemplateType> TemplateIteratorType;
	
		typedef  itk::ImageRegionConstIteratorWithIndex<ITK_TYPENAME Superclass::GradientImageType> GradientIteratorType;
	
	
		TemplateIteratorType ti( _template, this->GetFixedImageRegion() );
	
		typename FixedImageType::IndexType index;
	
		this->m_NumberOfPixelsCounted = 0;
	
		this->SetTransformParameters( parameters );
	
		const unsigned int ParametersDimension = this->GetNumberOfParameters();
		derivative = DerivativeType( ParametersDimension );
		derivative.Fill( itk::NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );
	
		ti.GoToBegin();
	
		for(ti.GoToBegin();!ti.IsAtEnd();++ti)
		{
      if(!ti.Value()) continue;
			index = ti.GetIndex();
			
			typename Superclass::InputPointType inputPoint;
			fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );
      
			typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );
      
      if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
			{
				const RealType movingValue = this->m_Interpolator->Evaluate(transformedPoint);
	
				const TransformJacobianType & jacobian = this->m_Transform->GetJacobian( inputPoint ); 
	
				
				const RealType fixedValue     =   fixedImage->GetPixel(index);
				this->m_NumberOfPixelsCounted++;
				const RealType diff = movingValue - fixedValue; 
	
				// Get the gradient by NearestNeighboorInterpolation: 
				// which is equivalent to round up the point components.
				typedef typename Superclass::OutputPointType OutputPointType;
				typedef typename OutputPointType::CoordRepType CoordRepType;
				typedef itk::ContinuousIndex< CoordRepType, MovingImageType::ImageDimension> MovingImageContinuousIndexType;
	
				MovingImageContinuousIndexType tempIndex;
				this->m_MovingImage->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );
	
				typename MovingImageType::IndexType mappedIndex; 
				for( unsigned int j = 0; j < MovingImageType::ImageDimension; j++ )
				{
					mappedIndex[j] = static_cast<long>( vnl_math_rnd( tempIndex[j] ) );
				}
	
				const GradientPixelType gradient = this->GetGradientImage()->GetPixel( mappedIndex );
	
				for(unsigned int par=0; par<ParametersDimension; par++)
				{
					RealType sum = itk::NumericTraits< RealType >::Zero;
					for(unsigned int dim=0; dim<ImageDimension; dim++)
					{
						sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
					}
					derivative[par] += sum;
				}
			}
		}
	
		if( !this->m_NumberOfPixelsCounted )
		{
			itkExceptionMacro(<<"All the points mapped to outside of the moving image");
		}
		else
		{
			for(unsigned int i=0; i<ParametersDimension; i++)
			{
				derivative[i] /= this->m_NumberOfPixelsCounted;
			}
		}
	}
	
	
	/*
	 * Get both the match Measure and theDerivative Measure 
	 */
	template <class TFixedImage, class TMovingImage, class TFixedTemplate> 
	void
	MeanSquaresImageToImageMetric<TFixedImage,TMovingImage,TFixedTemplate>
	::GetValueAndDerivative(const TransformParametersType & parameters, 
													MeasureType & value, DerivativeType  & derivative) const
	{
	
		itkExceptionMacro( << "Not Implemented" );
    //std::cout<<"*"<<std::flush;
		itkDebugMacro("GetValueAndDerivative( " << parameters << " ) ");
	
		if( !this->GetGradientImage() )
		{
			itkExceptionMacro(<<"The gradient image is null, maybe you forgot to call Initialize()");
		}
	
		FixedImageConstPointer fixedImage = this->m_FixedImage;
	
		if( !fixedImage ) 
		{
			itkExceptionMacro( << "Fixed image has not been assigned" );
		}
	
		const unsigned int ImageDimension = FixedImageType::ImageDimension;
	
		typedef  itk::ImageRegionConstIteratorWithIndex<
			FixedTemplateType> TemplateIteratorType;
	
		typedef  itk::ImageRegionConstIteratorWithIndex<
			ITK_TYPENAME Superclass::GradientImageType> GradientIteratorType;
	
	
		TemplateIteratorType ti( _template, this->GetFixedImageRegion() );
	
		typename FixedImageType::IndexType index;
	
		MeasureType measure = itk::NumericTraits< MeasureType >::Zero;
	
		this->m_NumberOfPixelsCounted = 0;
	
		this->SetTransformParameters( parameters );
	
		const unsigned int ParametersDimension = this->GetNumberOfParameters();
		derivative = DerivativeType( ParametersDimension );
		derivative.Fill( itk::NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );
	

		for(ti.GoToBegin();!ti.IsAtEnd();++ti)
		{
      if(!ti.Value()) continue;
      
			index = ti.GetIndex();
			
			typename Superclass::InputPointType inputPoint;
			fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );
	
			typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );
	
	
			if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
			{
				const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );
	
				const TransformJacobianType & jacobian =
					this->m_Transform->GetJacobian( inputPoint ); 
	
				
				const RealType fixedValue     = fixedImage->GetPixel(index);
				this->m_NumberOfPixelsCounted++;
	
				const RealType diff = movingValue - fixedValue; 
		
				measure += diff * diff;
	
				// Get the gradient by NearestNeighboorInterpolation: 
				// which is equivalent to round up the point components.
				typedef typename Superclass::OutputPointType OutputPointType;
				typedef typename OutputPointType::CoordRepType CoordRepType;
				typedef itk::ContinuousIndex< CoordRepType, MovingImageType::ImageDimension> 					MovingImageContinuousIndexType;
	
				MovingImageContinuousIndexType tempIndex;
				this->m_MovingImage->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );
	
				typename MovingImageType::IndexType mappedIndex; 
				for( unsigned int j = 0; j < MovingImageType::ImageDimension; j++ )
				{
					mappedIndex[j] = static_cast<long>( vnl_math_rnd( tempIndex[j] ) );
				}
	
				const GradientPixelType gradient = 
					this->GetGradientImage()->GetPixel( mappedIndex );
	
				for(unsigned int par=0; par<ParametersDimension; par++)
				{
					RealType sum = itk::NumericTraits< RealType >::Zero;
					for(unsigned int dim=0; dim<ImageDimension; dim++)
					{
						sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
					}
					derivative[par] += sum;
				}
			}
	
			++ti;
		}
	
		if( !this->m_NumberOfPixelsCounted )
		{
			itkExceptionMacro(<<"All the points mapped to outside of the moving image");
		}
		else
		{
			for(unsigned int i=0; i<ParametersDimension; i++)
			{
				derivative[i] /= this->m_NumberOfPixelsCounted;
			}
			measure /= this->m_NumberOfPixelsCounted;
		}
	
		value = measure;
	
	}

}; // end namespace minc
#endif
