/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanSquaresImageToImageMetric.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:35 $
  Version:   $Revision: 1.27 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMeanSquaresImageToImageMetric_h
#define __itkMeanSquaresImageToImageMetric_h

#include "itkImageToImageMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include <gsl/gsl_rng.h>

namespace minc
{
  //! Mean sqauares image difference metric, V.F. implementation
	template < class TFixedImage, class TMovingImage,class TFixedTemplate > 
	class MeanSquaresImageToImageMetric: 
			public itk::ImageToImageMetric< TFixedImage, TMovingImage>
	{
	public:
		/** Standard class typedefs. */
		typedef MeanSquaresImageToImageMetric    Self;
		typedef itk::ImageToImageMetric<TFixedImage, TMovingImage >  Superclass;
	
		typedef itk::SmartPointer<Self>         Pointer;
		typedef itk::SmartPointer<const Self>   ConstPointer;
    typedef TFixedTemplate                  FixedTemplateType;
	
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
	 
		/** Run-time type information (and related methods). */
		itkTypeMacro(MeanSquaresImageToImageMetric, itk::ImageToImageMetric);
	 
		/** Types transferred from the base class */
		typedef typename Superclass::RealType                 RealType;
		typedef typename Superclass::TransformType            TransformType;
		typedef typename Superclass::TransformPointer         TransformPointer;
		typedef typename Superclass::TransformParametersType  TransformParametersType;
		typedef typename Superclass::TransformJacobianType    TransformJacobianType;
		typedef typename Superclass::GradientPixelType        GradientPixelType;
	
		typedef typename Superclass::MeasureType              MeasureType;
		typedef typename Superclass::DerivativeType           DerivativeType;
		typedef typename Superclass::FixedImageType           FixedImageType;
		typedef typename Superclass::MovingImageType          MovingImageType;
		typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
		typedef typename Superclass::MovingImageConstPointer  MovingImageConstPointer;
		typedef typename TFixedTemplate::Pointer              FixedTemplatePointer;
		typedef typename TFixedTemplate::ConstPointer            FixedTemplateConstPointer;
    typedef typename Superclass::InterpolatorPointer  InterpolatorPointer;
    typedef typename Superclass::InterpolatorType     InterpolatorType;
		
    //! set template image
		virtual void SetTemplate(const TFixedTemplate * tmp) 
		{
			_template=tmp;
		}
    
    //! attach secondary images (not implemented)
		virtual void AttachSecondary(const FixedImageType* fixed,const TFixedTemplate * tmp, InterpolatorType * interpolator) 
		{
			_template2=tmp;
      _fixed2=fixed;
      _interpolator2=interpolator;
		}
	
		/** Get the derivatives of the match measure. */
		void GetDerivative( const TransformParametersType & parameters,
												DerivativeType & derivative ) const;
	
		/**  Get the value for single valued optimizers. */
		MeasureType GetValue( const TransformParametersType & parameters ) const;
    
    
		/**  Get value and derivatives for multiple valued optimizers. */
		void GetValueAndDerivative( const TransformParametersType & parameters,
																MeasureType& Value, DerivativeType& Derivative ) const;
	
    //! randomly subsample images, set ratio                                 
    void SetKeepRatio(double k) 
    {
      if(k<0.0) k=1.0;
      if(k>1.0) k=1.0;
      _keep=k;
    }
    
    double GetKeepRatio(void)
    {
      return _keep;
    }
	protected:
		MeanSquaresImageToImageMetric();
		virtual ~MeanSquaresImageToImageMetric();
		FixedTemplateConstPointer _template;
    
    //secondary targets
    FixedTemplateConstPointer _template2;
    FixedImageConstPointer _fixed2;
    InterpolatorPointer _interpolator2;
    
    //! use random subsampling if _keep<1.0 
    double _keep;
    gsl_rng * _rng;
    
    
    //! calculate mean square difference
		MeasureType GetMean( const TransformParametersType & parameters, const FixedImageType* fixed, const InterpolatorType* interpolator, const  FixedTemplateType* mask ) const;
	
	private:
		MeanSquaresImageToImageMetric(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
	
	};

} // end namespace itk

#include "mincMeanSquaresImageToImageMetric.txx"

#endif



