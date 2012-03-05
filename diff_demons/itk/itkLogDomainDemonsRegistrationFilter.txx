#ifndef __itkLogDomainDemonsRegistrationFilter_txx
#define __itkLogDomainDemonsRegistrationFilter_txx

#include "itkLogDomainDemonsRegistrationFilter.h"


namespace itk {

// Default constructor
template <class TFixedImage, class TMovingImage, class TField,class TMask>
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::LogDomainDemonsRegistrationFilterM()
{
    DemonsRegistrationFunctionPointer drfp = DemonsRegistrationFunctionType::New();
    this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
                                     drfp.GetPointer() ) );

    m_Multiplier = MultiplyByConstantType::New();
    m_Multiplier->InPlaceOn();

    m_BCHFilter = BCHFilterType::New();
    m_BCHFilter->InPlaceOn();

    // Set number of terms in the BCH approximation to default value
    m_BCHFilter->SetNumberOfApproximationTerms( 2 );
}


// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField,class TMask>
typename LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::DemonsRegistrationFunctionType*
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::DownCastDifferenceFunctionType()
{
    DemonsRegistrationFunctionType *drfp =
        dynamic_cast<DemonsRegistrationFunctionType *>(this->GetDifferenceFunction().GetPointer());

    if ( !drfp )
    {
        itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

    return drfp;
}


// Checks whether the DifferenceFunction is of type DemonsRegistrationFunction.
template <class TFixedImage, class TMovingImage, class TField,class TMask>
const typename LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::DemonsRegistrationFunctionType*
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::DownCastDifferenceFunctionType() const
{
    const DemonsRegistrationFunctionType *drfp =
        dynamic_cast<const DemonsRegistrationFunctionType *>(this->GetDifferenceFunction().GetPointer());

    if ( !drfp )
    {
        itkExceptionMacro( << "Could not cast difference function to SymmetricDemonsRegistrationFunction" );
    }

    return drfp;
}


// Set the function state values before each iteration
template <class TFixedImage, class TMovingImage, class TField,class TMask>
void
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::InitializeIteration()
{
    //std::cout<<"LogDomainDemonsRegistrationFilterM::InitializeIteration"<<std::endl;
    // update variables in the equation object
    DemonsRegistrationFunctionType *f = this->DownCastDifferenceFunctionType();
    f->SetDeformationField( this->GetDeformationField() );

    // call the superclass  implementation ( initializes f )
    Superclass::InitializeIteration();
}


// Get the metric value from the difference function
template <class TFixedImage, class TMovingImage, class TField,class TMask>
double
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::GetMetric() const
{
    const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    return drfp->GetMetric();
}


// Get Intensity Difference Threshold
template <class TFixedImage, class TMovingImage, class TField,class TMask>
double
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::GetIntensityDifferenceThreshold() const
{
    const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    return drfp->GetIntensityDifferenceThreshold();
}

// Set Intensity Difference Threshold
template <class TFixedImage, class TMovingImage, class TField,class TMask>
void
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::SetIntensityDifferenceThreshold(double threshold)
{
    DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    drfp->SetIntensityDifferenceThreshold(threshold);
}


// Set Maximum Update Step Length
template <class TFixedImage, class TMovingImage, class TField,class TMask>
void
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::SetMaximumUpdateStepLength(double step)
{
    DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    drfp->SetMaximumUpdateStepLength(step);
}

// Get Maximum Update Step Length
template <class TFixedImage, class TMovingImage, class TField,class TMask>
double
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::GetMaximumUpdateStepLength() const
{
    const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    return drfp->GetMaximumUpdateStepLength();
}


// Set number of terms used in the BCH approximation
template <class TFixedImage, class TMovingImage, class TField,class TMask>
void
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::SetNumberOfBCHApproximationTerms(unsigned int numterms)
{
    this->m_BCHFilter->SetNumberOfApproximationTerms(numterms);
}


// Get number of terms used in the BCH approximation
template <class TFixedImage, class TMovingImage, class TField,class TMask>
unsigned int
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::GetNumberOfBCHApproximationTerms() const
{
    return this->m_BCHFilter->GetNumberOfApproximationTerms();
}


// Get the metric value from the difference function
template <class TFixedImage, class TMovingImage, class TField,class TMask>
const double &
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::GetRMSChange() const
{
    const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    return drfp->GetRMSChange();
}


// Get gradient type
template <class TFixedImage, class TMovingImage, class TField,class TMask>
typename LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>::GradientType
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::GetUseGradientType() const
{
    const DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    return drfp->GetUseGradientType();
}

// Set gradient type
template <class TFixedImage, class TMovingImage, class TField,class TMask>
void
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::SetUseGradientType(GradientType gtype)
{
    DemonsRegistrationFunctionType *drfp = this->DownCastDifferenceFunctionType();
    drfp->SetUseGradientType(gtype);
}

// Get the metric value from the difference function
template <class TFixedImage, class TMovingImage, class TField,class TMask>
void
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::ApplyUpdate(
#if ( ITK_VERSION_MAJOR > 3 ) 
 const TimeStepType &dt
#else  
  TimeStepType dt
#endif  
)
{
  
    //std::cout<<"LogDomainDemonsRegistrationFilterM::ApplyUpdate"<<std::endl;
    // If we smooth the update buffer before applying it, then the are
    // approximating a viscuous problem as opposed to an elastic problem
    if ( this->GetSmoothUpdateField() )
    {
        this->SmoothUpdateField();
    }

    // Use time step if necessary. In many cases
    // the time step is one so this will be skipped
    if ( fabs(dt - 1.0)>1.0e-4 )
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


    // Apply update by using BCH approximation
    m_BCHFilter->SetInput( 0, this->GetOutput() );
    m_BCHFilter->SetInput( 1, this->GetUpdateBuffer() );
    if ( m_BCHFilter->GetInPlace() )
    {
        m_BCHFilter->GraftOutput( this->GetOutput() );
    }
    else
    {
        // Work-around for http://www.itk.org/Bug/view.php?id=8672
        m_BCHFilter->GraftOutput( DeformationFieldType::New() );
    }
    m_BCHFilter->GetOutput()->SetRequestedRegion( this->GetOutput()->GetRequestedRegion() );

    // Triggers in place update
    m_BCHFilter->Update();

    // Region passing stuff
    this->GraftOutput( m_BCHFilter->GetOutput() );


    //Smooth the velocity field
    if ( this->GetSmoothVelocityField() )
    {
        this->SmoothVelocityField();
    }

}


template <class TFixedImage, class TMovingImage, class TField,class TMask>
void
LogDomainDemonsRegistrationFilterM<TFixedImage,TMovingImage,TField,TMask>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf( os, indent );

    os << indent << "Multiplier: " << m_Multiplier << std::endl;
    os << indent << "BCHFilter: " << m_BCHFilter << std::endl;
}


} // end namespace itk

#endif
