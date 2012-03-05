/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsRegistrationFilterM.txx,v $
  Language:  C++
  Date:      $Date: 2008-12-08 16:00:52 $
  Version:   $Revision: 1.14 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsRegistrationFilterM_txx
#define __itkDemonsRegistrationFilterM_txx
#include "itkDemonsRegistrationFilterM.h"

namespace itk
{

/**
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
DemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::DemonsRegistrationFilterM()
{
 
  typename DemonsRegistrationFunctionType::Pointer drfp;
  drfp = DemonsRegistrationFunctionType::New();

  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                 drfp.GetPointer() ) );

  m_UseMovingImageGradient = false;

}


template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::PrintSelf(std::ostream& os, Indent indent) const
{ 
  Superclass::PrintSelf( os, indent );
  os << indent << "UseMovingImageGradient: ";
  os << m_UseMovingImageGradient << std::endl;
  os << indent << "Intensity difference threshold: " <<
    this->GetIntensityDifferenceThreshold() << std::endl;
}


/*
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::InitializeIteration()
{
 
  // call the superclass  implementation
  Superclass::InitializeIteration();

  // set the gradient selection flag
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << 
                       "Could not cast difference function to DemonsRegistrationFunction" );
    }
   
  drfp->SetUseMovingImageGradient( m_UseMovingImageGradient );


  /**
   * Smooth the deformation field
   */
  if ( this->GetSmoothDeformationField() )
    {
    this->SmoothDeformationField();
    }

}


/**
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
DemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetMetric() const
{
 
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << 
                       "Could not cast difference function to DemonsRegistrationFunction" );
    }
   
  return drfp->GetMetric();

}


/**
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
DemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetIntensityDifferenceThreshold() const
{
 
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << 
                       "Could not cast difference function to DemonsRegistrationFunction" );
    }
   
  return drfp->GetIntensityDifferenceThreshold();

}

/**
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::SetIntensityDifferenceThreshold(double threshold) 
{
 
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << 
                       "Could not cast difference function to DemonsRegistrationFunction" );
    }
   
  drfp->SetIntensityDifferenceThreshold(threshold);

}

/**
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::ApplyUpdate(
#if ( ITK_VERSION_MAJOR > 3 ) 
 const TimeStepType &dt
#else  
  TimeStepType dt
#endif  
)
{
  // If we smooth the update buffer before applying it, then the are
  // approximating a viscuous problem as opposed to an elastic problem
  if ( this->GetSmoothUpdateField() )
    {
    this->SmoothUpdateField();
    }
  
  this->Superclass::ApplyUpdate(dt);

  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
    (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << 
                       "Could not cast difference function to DemonsRegistrationFunction" );
    }

  this->SetRMSChange( drfp->GetRMSChange() );
   
}

} // end namespace itk

#endif
