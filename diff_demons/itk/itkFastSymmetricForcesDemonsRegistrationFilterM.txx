/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFastSymmetricForcesDemonsRegistrationFilterM.txx,v $
  Language:  C++
  Date:      $Date: 2009-04-05 23:09:19 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFastSymmetricForcesDemonsRegistrationFilterM_txx
#define __itkFastSymmetricForcesDemonsRegistrationFilterM_txx

#include "itkFastSymmetricForcesDemonsRegistrationFilterM.h"

namespace itk {

/**
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::FastSymmetricForcesDemonsRegistrationFilterM()
{
  
#if (ITK_VERSION_MAJOR < 4)
     this->SetNumberOfRequiredInputs(2);
#else
    //HACK: This really should define the names of the required inputs.
    this->SetNumberOfIndexedInputs(2);
    // Primary input is optional in this filter
    this->RemoveRequiredInputName( "Primary" );
#endif
  
  
  typename DemonsRegistrationFunctionType::Pointer drfp;
  drfp = DemonsRegistrationFunctionType::New();

  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                 drfp.GetPointer() ) );

  m_Multiplier = MultiplyByConstantType::New();
  m_Multiplier->InPlaceOn();

  m_Adder = AdderType::New();
  m_Adder->InPlaceOn();
}


/*
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::InitializeIteration()
{
  // update variables in the equation object
  DemonsRegistrationFunctionType *f = this->DownCastDifferenceFunctionType();

  f->SetDeformationField( this->GetDeformationField() );

  // call the superclass  implementation ( initializes f )
  Superclass::InitializeIteration();
}


/*
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetMetric() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetMetric();
}

/**
 * Return intensity difference threshold
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetIntensityDifferenceThreshold() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetIntensityDifferenceThreshold();
}

/**
 * Sets the intensity difference threshold
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::SetIntensityDifferenceThreshold(double threshold) 
{
  DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  drfp->SetIntensityDifferenceThreshold(threshold);
}


/**
 * Get the maximum update step length
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetMaximumUpdateStepLength() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetMaximumUpdateStepLength();
}

/**
 * Set the maximum update step length
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::SetMaximumUpdateStepLength(double threshold) 
{
  DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  drfp->SetMaximumUpdateStepLength(threshold);
}


/**
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
const double &
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetRMSChange() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetRMSChange();
}


/**
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
typename FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GradientType
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetUseGradientType() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetUseGradientType();
}

/**
 * 
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::SetUseGradientType(GradientType gtype) 
{
  DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  drfp->SetUseGradientType(gtype);
}

/**
 * Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
 * It throws and exception, if it is not.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
typename FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>::DemonsRegistrationFunctionType *
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::DownCastDifferenceFunctionType()
{
  DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << 
      "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
}
 
/**
 * Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
 * It throws and exception, if it is not.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
const typename FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>::DemonsRegistrationFunctionType *
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::DownCastDifferenceFunctionType() const
{
  const DemonsRegistrationFunctionType *drfp = 
    dynamic_cast<const DemonsRegistrationFunctionType *>
      (this->GetDifferenceFunction().GetPointer());
 
  if( !drfp )
    {
    itkExceptionMacro( << 
      "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

  return drfp;
}
 
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::AllocateUpdateBuffer()
{
  // The update buffer looks just like the output.
  DeformationFieldPointer output = this->GetOutput();
  DeformationFieldPointer upbuf = this->GetUpdateBuffer();

  upbuf->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  upbuf->SetRequestedRegion(output->GetRequestedRegion());
  upbuf->SetBufferedRegion(output->GetBufferedRegion());
  upbuf->SetOrigin(output->GetOrigin());
  upbuf->SetSpacing(output->GetSpacing());
  upbuf->SetDirection(output->GetDirection());
  upbuf->Allocate();
}


/**
 * Get the metric value from the difference function
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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

  // use time step if necessary
  if ( vcl_fabs(dt - 1.0)>1.0e-4 )
    {
    itkDebugMacro( "Using timestep: " << dt );
    m_Multiplier->SetConstant( dt );
    m_Multiplier->SetInput( this->GetUpdateBuffer() );
    m_Multiplier->GraftOutput( this->GetUpdateBuffer() );
    // in place update
    m_Multiplier->Update();
    // graft output back to this->GetUpdateBuffer()
    this->GetUpdateBuffer()->Graft( m_Multiplier->GetOutput() );
    }
  
  m_Adder->SetInput1( this->GetOutput() );
  m_Adder->SetInput2( this->GetUpdateBuffer() );

  m_Adder->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );
  m_Adder->Update();
  
  // Region passing stuff
  this->GraftOutput( m_Adder->GetOutput() );

  DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  
  this->SetRMSChange( drfp->GetRMSChange() );

  /*
   * Smooth the deformation field
   */
  if ( this->GetSmoothDeformationField() )
  {
    this->SmoothDeformationField();
  }
}

template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
FastSymmetricForcesDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::PrintSelf(std::ostream& os, Indent indent) const
{ 
  Superclass::PrintSelf( os, indent );
  os << indent << "Intensity difference threshold: " <<
    this->GetIntensityDifferenceThreshold() << std::endl;
}


} // end namespace itk

#endif
