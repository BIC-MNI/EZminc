/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFastSymmetricForcesDemonsRegistrationFilterM.h,v $
  Language:  C++
  Date:      $Date: 2008-12-06 13:28:10 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFastSymmetricForcesDemonsRegistrationFilterM_h
#define __itkFastSymmetricForcesDemonsRegistrationFilterM_h

#include "itkPDEDeformableRegistrationFilterM.h"
#include "itkESMDemonsRegistrationFunctionM.h"


#if ( ITK_VERSION_MAJOR > 3 ) 
#include <itkExponentialDisplacementFieldImageFilter.h>
#include <itkMultiplyImageFilter.h>
#else
#include <itkExponentialDeformationFieldImageFilter.h>
#include "itkMultiplyByConstantImageFilter.h"
#endif

#include "itkWarpVectorImageFilter.h"
#include "itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h"
#include "itkAddImageFilter.h"

namespace itk {

/** \class FastSymmetricForcesDemonsRegistrationFilterM
 * \brief Deformably register two images using a symmetric forces demons algorithm.
 *
 * This class was contributed by Tom Vercauteren, INRIA & Mauna Kea Technologies
 * based on a variation of the DemonsRegistrationFilter.
 *
 * FastSymmetricForcesDemonsRegistrationFilterM implements the demons deformable algorithm that 
 * register two images by computing the deformation field which will map a 
 * moving image onto a fixed image.
 *
 * A deformation field is represented as a image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the fixed image. The vector type must support element access via operator
 * []. It is assumed that the vector elements behave like floating point
 * scalars.
 *
 * This class is templated over the fixed image type, moving image type
 * and the deformation field type.
 *
 * The input fixed and moving images are set via methods SetFixedImage
 * and SetMovingImage respectively. An initial deformation field maybe set via
 * SetInitialDeformationField or SetInput. If no initial field is set,
 * a zero field is used as the initial condition.
 *
 * The output deformation field can be obtained via methods GetOutput
 * or GetDeformationField.
 *
 * This class make use of the finite difference solver hierarchy. Update
 * for each iteration is computed in DemonsRegistrationFunction.
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 * 
 * \sa DemonsRegistrationFilter
 * \sa DemonsRegistrationFunction
 * \ingroup DeformableImageRegistration MultiThreaded
 */
template<class TFixedImage, class TMovingImage, class TDeformationField,class TMask=TFixedImage>
class ITK_EXPORT FastSymmetricForcesDemonsRegistrationFilterM : 
    public PDEDeformableRegistrationFilterM< TFixedImage, TMovingImage,
                                            TDeformationField,TMask>
{
public:
  /** Standard class typedefs. */
  typedef FastSymmetricForcesDemonsRegistrationFilterM     Self;
  typedef PDEDeformableRegistrationFilterM<
    TFixedImage, TMovingImage,TDeformationField,TMask>          Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( FastSymmetricForcesDemonsRegistrationFilterM, 
                PDEDeformableRegistrationFilter );

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType             FixedImageType;
  typedef typename Superclass::FixedImagePointer          FixedImagePointer;

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType            MovingImageType;
  typedef typename Superclass::MovingImagePointer         MovingImagePointer;
  
  /** Mask image type. */
  typedef typename Superclass::MaskImageType            MaskImageType;
  typedef typename Superclass::MaskImagePointer         MaskImagePointer;  
  
  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType       DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer    DeformationFieldPointer;

  /** Get the metric value. The metric value is the mean square difference 
   * in intensity between the fixed image and transforming moving image 
   * computed over the the overlapping region between the two images. 
   * This value is calculated for the current iteration */
  virtual double GetMetric() const;
  virtual const double &GetRMSChange() const;

  /** DemonsRegistrationFilterFunction type. 
   * 
   *  FIXME: Why is this the only permissible function ?
   *
   */
  typedef ESMDemonsRegistrationFunctionM<
    FixedImageType, 
    MovingImageType, DeformationFieldType,MaskImageType>                DemonsRegistrationFunctionType;

  typedef typename DemonsRegistrationFunctionType::GradientType GradientType;
  virtual void SetUseGradientType( GradientType gtype );
  virtual GradientType GetUseGradientType() const;

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void SetIntensityDifferenceThreshold(double);
  virtual double GetIntensityDifferenceThreshold() const;

  virtual void SetMaximumUpdateStepLength(double);
  virtual double GetMaximumUpdateStepLength() const;

protected:
  FastSymmetricForcesDemonsRegistrationFilterM();
  ~FastSymmetricForcesDemonsRegistrationFilterM() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize the state of filter and equation before each iteration. */
  virtual void InitializeIteration();

  /** This method allocates storage in m_UpdateBuffer.  It is called from
   * FiniteDifferenceFilter::GenerateData(). */
  virtual void AllocateUpdateBuffer();

  /** FiniteDifferenceFunction type. */
  typedef typename 
    Superclass::FiniteDifferenceFunctionType              FiniteDifferenceFunctionType;

  /** Take timestep type from the FiniteDifferenceFunction. */
  typedef typename 
    FiniteDifferenceFunctionType::TimeStepType            TimeStepType;

  /** Apply update. */
#if ( ITK_VERSION_MAJOR > 3 ) 
  virtual void ApplyUpdate(const TimeStepType &dt);  
  typedef MultiplyImageFilter<
    DeformationFieldType, 
    DeformationFieldType, DeformationFieldType >         MultiplyByConstantType;

#else  
  virtual void ApplyUpdate(TimeStepType dt);
  /** other typedefs */
  typedef MultiplyByConstantImageFilter<
    DeformationFieldType, 
    TimeStepType, DeformationFieldType >                  MultiplyByConstantType;
#endif  


  typedef AddImageFilter<
     DeformationFieldType, 
     DeformationFieldType, DeformationFieldType>          AdderType;


  typedef typename MultiplyByConstantType::Pointer        MultiplyByConstantPointer;
  typedef typename AdderType::Pointer                     AdderPointer;


private:
  FastSymmetricForcesDemonsRegistrationFilterM(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Downcast the DifferenceFunction using a dynamic_cast to ensure that it is of the correct type.
   * this method will throw an exception if the function is not of the expected type. */
  DemonsRegistrationFunctionType *  DownCastDifferenceFunctionType();
  const DemonsRegistrationFunctionType *  DownCastDifferenceFunctionType() const;

  MultiplyByConstantPointer         m_Multiplier;
  AdderPointer                      m_Adder;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFastSymmetricForcesDemonsRegistrationFilterM.txx"
#endif

#endif
