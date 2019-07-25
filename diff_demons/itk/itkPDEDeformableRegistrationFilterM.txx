/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPDEDeformableRegistrationFilterM.txx,v $
  Language:  C++
  Date:      $Date: 2009-01-26 21:45:56 $
  Version:   $Revision: 1.32 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPDEDeformableRegistrationFilterM_txx
#define __itkPDEDeformableRegistrationFilterM_txx

#include "itkPDEDeformableRegistrationFilterM.h"

#include "itkExceptionObject.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkDataObject.h"

#include "itkGaussianOperator.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"

#include "itkMath.h"

namespace itk
  {

  /**
   * Default constructor
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::PDEDeformableRegistrationFilterM()
  {

    this->SetNumberOfRequiredInputs ( 2 );

    this->SetNumberOfIterations ( 10 );

    unsigned int j;
    for ( j = 0; j < ImageDimension; j++ )
      {
        m_StandardDeviations[j] = 1.0;
        m_UpdateFieldStandardDeviations[j] = 1.0;
      }

    m_TempField = DeformationFieldType::New();
    m_MaximumError = 0.1;
    m_MaximumKernelWidth = 30;
    m_StopRegistrationFlag = false;

    m_SmoothDeformationField = true;
    m_SmoothUpdateField = false;

    m_MovingImageMask = NULL;
    m_FixedImageMask = NULL;
  }


  /*
   * Set the fixed image.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::SetFixedImage (
    const FixedImageType * ptr )
  {
    this->ProcessObject::SetNthInput ( 1, const_cast< FixedImageType * > ( ptr ) );
  }


  /*
   * Get the fixed image.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  const typename PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::FixedImageType *
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::GetFixedImage() const
    {
      return dynamic_cast< const FixedImageType * >
             ( this->ProcessObject::GetInput ( 1 ) );
    }


  /*
   * Set the moving image.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::SetMovingImage (
    const MovingImageType * ptr )
  {
    this->ProcessObject::SetNthInput ( 2, const_cast< MovingImageType * > ( ptr ) );
  }


  /*
   * Get the moving image.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  const typename PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::MovingImageType *
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::GetMovingImage() const
    {
      return dynamic_cast< const MovingImageType * >
             ( this->ProcessObject::GetInput ( 2 ) );
    }


  /*
   *
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  std::vector<SmartPointer<DataObject> >::size_type
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::GetNumberOfValidRequiredInputs() const
    {
      typename std::vector<SmartPointer<DataObject> >::size_type num = 0;

      if ( this->GetFixedImage() )
        {
          num++;
        }

      if ( this->GetMovingImage() )
        {
          num++;
        }

      return num;
    }


  /**
   * Set the standard deviations.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::SetStandardDeviations (
    double value )
  {

    unsigned int j;
    for ( j = 0; j < ImageDimension; j++ )
      {
        if ( value != m_StandardDeviations[j] )
          {
            break;
          }
      }
    if ( j < ImageDimension )
      {
        this->Modified();
        for ( j = 0; j < ImageDimension; j++ )
          {
            m_StandardDeviations[j] = value;
          }
      }

  }

  /*
   * Set the standard deviations.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::SetUpdateFieldStandardDeviations (
    double value )
  {

    unsigned int j;
    for ( j = 0; j < ImageDimension; j++ )
      {
        if ( value != m_UpdateFieldStandardDeviations[j] )
          {
            break;
          }
      }
    if ( j < ImageDimension )
      {
        this->Modified();
        for ( j = 0; j < ImageDimension; j++ )
          {
            m_UpdateFieldStandardDeviations[j] = value;
          }
      }

  }


  /*
   * Standard PrintSelf method.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::PrintSelf ( std::ostream& os, Indent indent ) const
    {
      Superclass::PrintSelf ( os, indent );
      
      os << indent << "MovingImageMask: ";
      os << m_MovingImageMask.GetPointer() << std::endl;
      os << indent << "FixedImageMask: ";
      os << m_FixedImageMask.GetPointer() << std::endl;
      
      os << indent << "Smooth deformation field: "
      << ( m_SmoothDeformationField ? "on" : "off" ) << std::endl;
      os << indent << "Standard deviations: [";
      unsigned int j;
      for ( j = 0; j < ImageDimension - 1; j++ )
        {
          os << m_StandardDeviations[j] << ", ";
        }
      os << m_StandardDeviations[j] << "]" << std::endl;
      os << indent << "Smooth update field: "
      << ( m_SmoothUpdateField ? "on" : "off" ) << std::endl;
      os << indent << "Update field standard deviations: [";
      for ( j = 0; j < ImageDimension - 1; j++ )
        {
          os << m_UpdateFieldStandardDeviations[j] << ", ";
        }
      os << m_UpdateFieldStandardDeviations[j] << "]" << std::endl;
      os << indent << "StopRegistrationFlag: ";
      os << m_StopRegistrationFlag << std::endl;
      os << indent << "MaximumError: ";
      os << m_MaximumError << std::endl;
      os << indent << "MaximumKernelWidth: ";
      os << m_MaximumKernelWidth << std::endl;
      
    }


  /*
   * Set the function state values before each iteration
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::InitializeIteration()
  {

    MovingImageConstPointer movingPtr = this->GetMovingImage();
    FixedImageConstPointer fixedPtr = this->GetFixedImage();

    if ( !movingPtr || !fixedPtr )
      {
        itkExceptionMacro ( << "Fixed and/or moving image not set" );
      }

    // update variables in the equation object
    PDEDeformableRegistrationFunctionType *f =
      dynamic_cast<PDEDeformableRegistrationFunctionType *>
      ( this->GetDifferenceFunction().GetPointer() );

    if ( !f )
      {
        itkExceptionMacro ( <<"FiniteDifferenceFunction not of type PDEDeformableRegistrationFilterMFunction" );
      }

    f->SetFixedImage ( fixedPtr );
    f->SetMovingImage ( movingPtr );
    
    if(m_MovingImageMask)
      f->SetMovingImageMask(m_MovingImageMask);
    
    if(m_FixedImageMask)
      f->SetFixedImageMask(m_FixedImageMask);
    

    this->Superclass::InitializeIteration();

  }


  /*
   * Override the default implemenation for the case when the
   * initial deformation is not set.
   * If the initial deformation is not set, the output is
   * fill with zero vectors.
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::CopyInputToOutput()
  {

    typename Superclass::InputImageType::ConstPointer  inputPtr  = this->GetInput();

    if ( inputPtr )
      {
        this->Superclass::CopyInputToOutput();
      }
    else
      {
        typename Superclass::PixelType zeros;
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
            zeros[j] = 0;
          }

        typename OutputImageType::Pointer output = this->GetOutput();

        ImageRegionIterator<OutputImageType> out ( output, output->GetRequestedRegion() );

        while ( ! out.IsAtEnd() )
          {
            out.Value() =  zeros;
            ++out;
          }
      }
  }


  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::GenerateOutputInformation()
  {

    typename DataObject::Pointer output;

    if ( this->GetInput ( 0 ) )
      {
        // Initial deformation field is set.
        // Copy information from initial field.
        this->Superclass::GenerateOutputInformation();

      }
    else if ( this->GetFixedImage() )
      {
        // Initial deforamtion field is not set.
        // Copy information from the fixed image.
        for ( unsigned int idx = 0; idx <
              this->GetNumberOfOutputs(); ++idx )
          {
            output = this->GetOutput ( idx );
            if ( output )
              {
                output->CopyInformation ( this->GetFixedImage() );
              }
          }

      }

  }


  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::GenerateInputRequestedRegion()
  {

    // call the superclass's implementation
    Superclass::GenerateInputRequestedRegion();

    // request the largest possible region for the moving image
    MovingImagePointer movingPtr =
      const_cast< MovingImageType * > ( this->GetMovingImage() );
    if ( movingPtr )
      {
        movingPtr->SetRequestedRegionToLargestPossibleRegion();
      }

    // just propagate up the output requested region for
    // the fixed image and initial deformation field.
    DeformationFieldPointer inputPtr =
      const_cast< DeformationFieldType * > ( this->GetInput() );
    DeformationFieldPointer outputPtr = this->GetOutput();
    FixedImagePointer fixedPtr =
      const_cast< FixedImageType *> ( this->GetFixedImage() );

    if ( inputPtr )
    {
        inputPtr->SetRequestedRegion ( outputPtr->GetRequestedRegion() );
    }

    if ( fixedPtr )
    {
        fixedPtr->SetRequestedRegion ( outputPtr->GetRequestedRegion() );
    }
    
    if( m_MovingImageMask )
    {
      MaskImageType* movingMaskPtr =
	const_cast< MaskImageType * > ( m_MovingImageMask.GetPointer() );
	
      movingMaskPtr->SetRequestedRegionToLargestPossibleRegion();
    }
    
    if( m_FixedImageMask )
    {
      MaskImageType* fixedMaskPtr =
	const_cast< MaskImageType * > ( m_FixedImageMask.GetPointer() );
	
      fixedMaskPtr->SetRequestedRegion ( outputPtr->GetRequestedRegion() );
    }

  }


  /*
   * Release memory of internal buffers
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::PostProcessOutput()
  {
    this->Superclass::PostProcessOutput();
    m_TempField->Initialize();
  }


  /*
   * Initialize flags
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::Initialize()
  {
    this->Superclass::Initialize();
    m_StopRegistrationFlag = false;
  }


  /*
   * Smooth deformation using a separable Gaussian kernel
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::SmoothDeformationField() //TODO: replace to 'regularize' deformation field
  {
    //std::cout<<"PDEDeformableRegistrationFilterM::SmoothDeformationField()"<<std::endl;
    DeformationFieldPointer field = this->GetOutput();

    // copy field to TempField
    m_TempField->SetOrigin ( field->GetOrigin() );
    m_TempField->SetSpacing ( field->GetSpacing() );
    m_TempField->SetDirection ( field->GetDirection() );
    m_TempField->SetLargestPossibleRegion (
      field->GetLargestPossibleRegion() );
    m_TempField->SetRequestedRegion (
      field->GetRequestedRegion() );
    m_TempField->SetBufferedRegion ( field->GetBufferedRegion() );
    m_TempField->Allocate();

    typedef typename DeformationFieldType::PixelType    VectorType;
    typedef typename VectorType::ValueType              ScalarType;
    typedef GaussianOperator<ScalarType,ImageDimension> OperatorType;
    typedef VectorNeighborhoodOperatorImageFilter<
    DeformationFieldType,
    DeformationFieldType>                             SmootherType;

    OperatorType * oper = new OperatorType;
    typename SmootherType::Pointer smoother = SmootherType::New();

    typedef typename DeformationFieldType::PixelContainerPointer
    PixelContainerPointer;
    PixelContainerPointer swapPtr;

    // graft the output field onto the mini-pipeline
    smoother->GraftOutput ( m_TempField );

    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
        // smooth along this dimension
        oper->SetDirection ( j );
        double variance = itk::Math::sqr ( m_StandardDeviations[j] );
        oper->SetVariance ( variance );
        oper->SetMaximumError ( m_MaximumError );
        oper->SetMaximumKernelWidth ( m_MaximumKernelWidth );
        oper->CreateDirectional();

        // todo: make sure we only smooth within the buffered region
        smoother->SetOperator ( *oper );
        smoother->SetInput ( field );
        smoother->Update();

        if ( j < ImageDimension - 1 )
          {
            // swap the containers
            swapPtr = smoother->GetOutput()->GetPixelContainer();
            smoother->GraftOutput ( field );
            field->SetPixelContainer ( swapPtr );
            smoother->Modified();
          }

      }

    // graft the output back to this filter
    m_TempField->SetPixelContainer ( field->GetPixelContainer() );
    this->GraftOutput ( smoother->GetOutput() );

    delete oper;

  }

  /*
   * Smooth deformation using a separable Gaussian kernel
   */
  template <class TFixedImage, class TMovingImage, class TDeformationField,class TMask>
  void
  PDEDeformableRegistrationFilterM<TFixedImage,TMovingImage,TDeformationField,TMask>
  ::SmoothUpdateField() //TODO: replace to 'regularize' update field
  {
    //std::cout<<"PDEDeformableRegistrationFilterM::SmoothUpdateField()"<<std::endl;
    // The update buffer will be overwritten with new data.
    DeformationFieldPointer field = this->GetUpdateBuffer();

    typedef typename DeformationFieldType::PixelType    VectorType;
    typedef typename VectorType::ValueType              ScalarType;
    typedef GaussianOperator<ScalarType,ImageDimension> OperatorType;
    typedef VectorNeighborhoodOperatorImageFilter<
      DeformationFieldType,
      DeformationFieldType>                             SmootherType;

    OperatorType opers[ImageDimension];
    
    typename SmootherType::Pointer smoothers[ImageDimension];

    for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
        // smooth along this dimension
        opers[j].SetDirection ( j );
        double variance = itk::Math::sqr ( this->GetUpdateFieldStandardDeviations() [j] );
        //double variance = itk::Math::sqr ( 1.0 );
        opers[j].SetVariance ( variance );
        opers[j].SetMaximumError ( this->GetMaximumError() );
        opers[j].SetMaximumKernelWidth ( this->GetMaximumKernelWidth() );
        opers[j].CreateDirectional();

        smoothers[j] = SmootherType::New();
        smoothers[j]->SetOperator ( opers[j] );
        smoothers[j]->ReleaseDataFlagOn();

        if ( j > 0 )
        {
            smoothers[j]->SetInput ( smoothers[j-1]->GetOutput() );
        }
    }
    
    smoothers[0]->SetInput ( field );
    smoothers[ImageDimension-1]->GetOutput() ->SetRequestedRegion ( field->GetBufferedRegion() );

    smoothers[ImageDimension-1]->Update();

    // field to contain the final smoothed data, do the equivalent of a graft
    field->SetPixelContainer ( smoothers[ImageDimension-1]->GetOutput()->GetPixelContainer() );
    field->SetRequestedRegion ( smoothers[ImageDimension-1]->GetOutput()->GetRequestedRegion() );
    field->SetBufferedRegion ( smoothers[ImageDimension-1]->GetOutput()->GetBufferedRegion() );
    field->SetLargestPossibleRegion ( smoothers[ImageDimension-1]->GetOutput()->GetLargestPossibleRegion() );
    
    field->CopyInformation ( smoothers[ImageDimension-1]->GetOutput() );
  }


} // end namespace itk

#endif
