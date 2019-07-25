/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffeomorphicDemonsRegistrationFilterM.txx,v $
  Language:  C++
  Date:      $Date: 2009-10-29 15:03:32 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDiffeomorphicDemonsRegistrationFilterM_txx
#define __itkDiffeomorphicDemonsRegistrationFilterM_txx

#include "itkDiffeomorphicDemonsRegistrationFilterM.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

namespace itk {

/**
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::DiffeomorphicDemonsRegistrationFilterM()
   :m_UseFirstOrderExp(false)
{
 
  typename DemonsRegistrationFunctionType::Pointer drfp;
  drfp = DemonsRegistrationFunctionType::New();

  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                 drfp.GetPointer() ) );

#if (ITK_VERSION_MAJOR < 4)
  this->SetNumberOfRequiredInputs(2);
#else
  //HACK: This really should define the names of the required inputs.
  this->SetNumberOfIndexedInputs(2);
  // Primary input is optional in this filter
  this->RemoveRequiredInputName( "Primary" );
#endif  
  
  m_Multiplier = MultiplyByConstantType::New();
  m_Multiplier->InPlaceOn();

  m_Exponentiator = FieldExponentiatorType::New();
  
  m_Warper = VectorWarperType::New();
  FieldInterpolatorPointer VectorInterpolator =
     FieldInterpolatorType::New();
  m_Warper->SetInterpolator(VectorInterpolator);

  m_Adder = AdderType::New();
  m_Adder->InPlaceOn();
}


/**
 * Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
 * It throws and exception, if it is not.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
typename DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>::DemonsRegistrationFunctionType *
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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
const typename DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>::DemonsRegistrationFunctionType *
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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
 

/**
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetMetric() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetMetric();
}

/**
 *  Get Intensity Difference Threshold
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetIntensityDifferenceThreshold() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetIntensityDifferenceThreshold();
}

/**
 *  Set Intensity Difference Threshold
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::SetIntensityDifferenceThreshold(double threshold) 
{
  DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  drfp->SetIntensityDifferenceThreshold(threshold);
}


/**
 *  Get Maximum Update Step Length
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
double
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetMaximumUpdateStepLength() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetMaximumUpdateStepLength();
}

/**
 *  Set Maximum Update Step Length
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GetRMSChange() const
{
  const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  return drfp->GetRMSChange();
}


/**
 *   
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
typename DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::GradientType
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::SetUseGradientType(GradientType gtype) 
{
  DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
  drfp->SetUseGradientType(gtype);
}


/**
 *
 */
template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
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

  // Use time step if necessary. In many cases
  // the time step is one so this will be skipped
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

  

  if ( this->m_UseFirstOrderExp )
    {
    // use s <- s o (Id +u)
    
    // skip exponential and compose the vector fields
    m_Warper->SetOutputOrigin( this->GetUpdateBuffer()->GetOrigin() );
    m_Warper->SetOutputSpacing( this->GetUpdateBuffer()->GetSpacing() );
    m_Warper->SetOutputDirection( this->GetUpdateBuffer()->GetDirection() );
    m_Warper->SetInput( this->GetOutput() );
  
    
    m_Adder->SetInput1( m_Warper->GetOutput() );
    m_Adder->SetInput2( this->GetUpdateBuffer() );
    
    m_Adder->GetOutput()->SetRequestedRegion(
       this->GetOutput()->GetRequestedRegion() );
    }
  else
    {
    // use s <- s o exp(u)
    
    // compute the exponential
    m_Exponentiator->SetInput( this->GetUpdateBuffer() );
    
    const double imposedMaxUpStep = this->GetMaximumUpdateStepLength();
    if( imposedMaxUpStep > 0.0 )
      {
      // max(norm(Phi))/2^N <= 0.25*pixelspacing
      const double numiterfloat = 2.0 + vcl_log(imposedMaxUpStep)/itk::Math::ln2;
      unsigned int numiter = 0;
      if ( numiterfloat > 0.0 )
        {
        numiter = Math::Ceil<unsigned int>( numiterfloat );
        }
      
      m_Exponentiator->AutomaticNumberOfIterationsOff();
      m_Exponentiator->SetMaximumNumberOfIterations( numiter );
      }
    else
      {
      m_Exponentiator->AutomaticNumberOfIterationsOn();
      // just set a high value so that automatic number of step
      // is not thresholded
      m_Exponentiator->SetMaximumNumberOfIterations( 2000u );
      }

    m_Exponentiator->GetOutput()->SetRequestedRegion(
      this->GetOutput()->GetRequestedRegion() );

    m_Exponentiator->Update();

    // compose the vector fields
    m_Warper->SetOutputOrigin( this->GetUpdateBuffer()->GetOrigin() );
    m_Warper->SetOutputSpacing( this->GetUpdateBuffer()->GetSpacing() );
    m_Warper->SetOutputDirection( this->GetUpdateBuffer()->GetDirection() );
    m_Warper->SetInput( this->GetOutput() );

#if ( ITK_VERSION_MAJOR > 3 ) 
    m_Warper->SetDisplacementField( m_Exponentiator->GetOutput() );
#else
    m_Warper->SetDeformationField( m_Exponentiator->GetOutput() );
#endif

    m_Warper->Update();
    
    m_Adder->SetInput1( m_Warper->GetOutput() );
    m_Adder->SetInput2( m_Exponentiator->GetOutput() );
    
    m_Adder->GetOutput()->SetRequestedRegion(
       this->GetOutput()->GetRequestedRegion() );
    }

  // Triggers update
  m_Adder->Update();

  // Region passing stuff
  this->GraftOutput( m_Adder->GetOutput() );
  
  DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();

  this->SetRMSChange( drfp->GetRMSChange() );

  /**
   * Smooth the deformation field
   */
  if( this->GetSmoothDeformationField() )
    {
    this->SmoothDeformationField();
    }
}

template <class TFixedImage, class TMovingImage, class TDeformationField, class TMask>
void
DiffeomorphicDemonsRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
::PrintSelf(std::ostream& os, Indent indent) const
{ 
  Superclass::PrintSelf( os, indent );
  os << indent << "Intensity difference threshold: " <<
    this->GetIntensityDifferenceThreshold() << std::endl;
  os << indent << "Use First Order exponential: " << 
    this->m_UseFirstOrderExp << std::endl;
}


} // end namespace itk

#endif
