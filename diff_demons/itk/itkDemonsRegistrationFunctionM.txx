/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsRegistrationFunctionM.txx,v $
  Language:  C++
  Date:      $Date: 2008-12-08 16:00:52 $
  Version:   $Revision: 1.32 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsRegistrationFunctionM_txx
#define __itkDemonsRegistrationFunctionM_txx

#include "itkDemonsRegistrationFunctionM.h"
#include "itkExceptionObject.h"
#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::DemonsRegistrationFunctionM()
{

  RadiusType r;
  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    r[j] = 0;
    }
  this->SetRadius(r);

  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  this->SetMovingImage(NULL);
  this->SetFixedImage(NULL);
  //m_FixedImageSpacing.Fill( 1.0 );
  //m_FixedImageOrigin.Fill( 0.0 );
  m_Normalizer = 1.0;
  m_FixedImageGradientCalculator = GradientCalculatorType::New();


  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );
  
  m_MovingImageMaskInterpolator = static_cast<InterpolatorType*>(DefaultMaskInterpolatorType::New());

  m_Metric = NumericTraits<double>::max();
  m_SumOfSquaredDifference = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_RMSChange = NumericTraits<double>::max();
  m_SumOfSquaredChange = 0.0;

  m_MovingImageGradientCalculator = MovingImageGradientCalculatorType::New();
  m_UseMovingImageGradient = false;

}


/**
 * Standard "PrintSelf" method.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "MovingImageIterpolator: ";
  os << m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "MovingImageMaskIterpolator: ";
  os << m_MovingImageMaskInterpolator.GetPointer() << std::endl;
  
  os << indent << "FixedImageGradientCalculator: ";
  os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
  os << indent << "DenominatorThreshold: ";
  os << m_DenominatorThreshold << std::endl;
  os << indent << "IntensityDifferenceThreshold: ";
  os << m_IntensityDifferenceThreshold << std::endl;

  os << indent << "UseMovingImageGradient: ";
  os << m_UseMovingImageGradient << std::endl;

  os << indent << "Metric: ";
  os << m_Metric << std::endl;
  os << indent << "SumOfSquaredDifference: ";
  os << m_SumOfSquaredDifference << std::endl;
  os << indent << "NumberOfPixelsProcessed: ";
  os << m_NumberOfPixelsProcessed << std::endl;
  os << indent << "RMSChange: ";
  os << m_RMSChange << std::endl;
  os << indent << "SumOfSquaredChange: ";
  os << m_SumOfSquaredChange << std::endl;

}


/**
 *
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::SetIntensityDifferenceThreshold(double threshold)
{
  m_IntensityDifferenceThreshold = threshold;
}

/**
 *
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetIntensityDifferenceThreshold() const
{
  return m_IntensityDifferenceThreshold;
}


/**
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::InitializeIteration()
{
  if( !this->GetMovingImage() || !this->GetFixedImage() || !m_MovingImageInterpolator )
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not set" );
    }

  // cache fixed image information
  SpacingType fixedImageSpacing    = this->GetFixedImage()->GetSpacing();
  m_ZeroUpdateReturn.Fill(0.0);

  // compute the normalizer
  m_Normalizer      = 0.0;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    m_Normalizer += fixedImageSpacing[k] * fixedImageSpacing[k];
    }
  m_Normalizer /= static_cast<double>( ImageDimension );


  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage( this->GetFixedImage() );
  m_MovingImageGradientCalculator->SetInputImage( this->GetMovingImage() );

  // setup moving image interpolator
  m_MovingImageInterpolator->SetInputImage( this->GetMovingImage() );

  
  if (this->GetMovingImageMask()) //we are using masks
  {
      m_MovingImageMaskInterpolator->SetInputImage( this->GetMovingImageMask() );

  }
  
  // initialize metric computation variables
  m_SumOfSquaredDifference  = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_SumOfSquaredChange      = 0.0;

}


/**
 * Compute update at a specify neighbourhood
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
typename DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::PixelType
DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::ComputeUpdate(const NeighborhoodType &it, void * gd,
                const FloatOffsetType& itkNotUsed(offset))
{
  // Get fixed image related information
  // Note: no need to check the index is within
  // fixed image buffer. This is done by the external filter.
  const IndexType index = it.GetIndex();
  const double fixedValue = (double) this->GetFixedImage()->GetPixel( index );

  // Get moving image related information
  PointType mappedPoint;
  this->GetFixedImage()->TransformIndexToPhysicalPoint(index, mappedPoint);
  
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    mappedPoint[j] += it.GetCenterPixel()[j];
    }

  MaskPixelType fixedMask=NumericTraits <MaskPixelType>::OneValue();
  if ( this->GetFixedImageMask() )
  {
      fixedMask = this->GetFixedImageMask()->GetPixel ( index );

      if ( fixedMask == NumericTraits <MaskPixelType>::ZeroValue() )
	  return m_ZeroUpdateReturn;
  }

  MaskPixelType movingMask=NumericTraits <MaskPixelType>::OneValue();

  if ( this->GetMovingImageMask() )
  {
    if(m_MovingImageMaskInterpolator->IsInsideBuffer( mappedPoint ) )
      movingMask = m_MovingImageMaskInterpolator->Evaluate ( mappedPoint );
  }

  if ( movingMask == NumericTraits <MaskPixelType>::ZeroValue() )
  {
      return m_ZeroUpdateReturn;
  }


  double movingValue;
  if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
    {
    movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
    }
  else
    {
    return m_ZeroUpdateReturn;
    }

  CovariantVectorType gradient;
  // Compute the gradient of either fixed or moving image
  if( !m_UseMovingImageGradient )
    {
    gradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
    }
  else
    {
    gradient = m_MovingImageGradientCalculator->Evaluate( mappedPoint );
    }

  double gradientSquaredMagnitude = 0;
  for(unsigned int j = 0; j < ImageDimension; j++ )
    {
    gradientSquaredMagnitude += vnl_math_sqr( gradient[j] );
    }

  /**
   * Compute Update.
   * In the original equation the denominator is defined as (g-f)^2 + grad_mag^2.
   * However there is a mismatch in units between the two terms.
   * The units for the second term is intensity^2/mm^2 while the
   * units for the first term is intensity^2. This mismatch is particularly
   * problematic when the fixed image does not have unit spacing.
   * In this implemenation, we normalize the first term by a factor K,
   * such that denominator = (g-f)^2/K + grad_mag^2
   * where K = mean square spacing to compensate for the mismatch in units.
   */
  const double speedValue = fixedValue - movingValue;
  const double sqr_speedValue=vnl_math_sqr(speedValue);

  // update the metric
  GlobalDataStruct *globalData = (GlobalDataStruct *)gd;
  if ( globalData )
    {
    globalData->m_SumOfSquaredDifference += sqr_speedValue;
    globalData->m_NumberOfPixelsProcessed += 1;
    }

  const double denominator = sqr_speedValue / m_Normalizer +
    gradientSquaredMagnitude;

  if ( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold ||
       denominator < m_DenominatorThreshold )
    {
    return m_ZeroUpdateReturn;
    }

  PixelType update;
  for(unsigned int j = 0; j < ImageDimension; j++ )
    {
    update[j] = speedValue * gradient[j] / denominator;
    if ( globalData )
      {
      globalData->m_SumOfSquaredChange += vnl_math_sqr( update[j] );
      }
    }
  return update;
}

/**
 * Update the metric and release the per-thread-global data.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFunctionM<TFixedImage,TMovingImage,TDeformationField,TMask>
::ReleaseGlobalDataPointer( void *gd ) const
{
  GlobalDataStruct * globalData = (GlobalDataStruct *) gd;

  m_MetricCalculationLock.Lock();
  m_SumOfSquaredDifference += globalData->m_SumOfSquaredDifference;
  m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
  m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
  if ( m_NumberOfPixelsProcessed )
    {
    m_Metric = m_SumOfSquaredDifference / 
      static_cast<double>( m_NumberOfPixelsProcessed ); 
    m_RMSChange = vcl_sqrt( m_SumOfSquaredChange / 
                            static_cast<double>( m_NumberOfPixelsProcessed ) ); 
    }
  m_MetricCalculationLock.Unlock();

  delete globalData;
}

} // end namespace itk

#endif
