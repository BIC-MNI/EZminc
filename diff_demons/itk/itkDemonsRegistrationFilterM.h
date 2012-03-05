/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDemonsRegistrationFilterM.h,v $
  Language:  C++
  Date:      $Date: 2009-04-23 03:53:35 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDemonsRegistrationFilterM_h
#define __itkDemonsRegistrationFilterM_h

#include "itkPDEDeformableRegistrationFilterM.h"
#include "itkDemonsRegistrationFunctionM.h"

namespace itk {

/** \class DemonsRegistrationFilterM
 * \brief Deformably register two images using the demons algorithm.
 *
 * DemonsRegistrationFilterM implements the demons deformable algorithm that 
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
 * The algorithm has one parameters: the number of iteration to be performed.
 *
 * The output deformation field can be obtained via methods GetOutput
 * or GetDeformationField.
 *
 * This class make use of the finite difference solver hierarchy. Update
 * for each iteration is computed in DemonsRegistrationFunction.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 * 
 * \sa DemonsRegistrationFunction 
 * \ingroup DeformableImageRegistration MultiThreaded
 */
template<class TFixedImage, class TMovingImage, class TDeformationField,class TMask=TFixedImage>
class ITK_EXPORT DemonsRegistrationFilterM : 
    public PDEDeformableRegistrationFilterM< TFixedImage, TMovingImage,
                                            TDeformationField, TMask>
{
public:
  /** Standard class typedefs. */
  typedef DemonsRegistrationFilterM                  Self;
  typedef PDEDeformableRegistrationFilterM<
    TFixedImage, TMovingImage,TDeformationField,TMask>    Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( DemonsRegistrationFilterM, 
                PDEDeformableRegistrationFilterM );

  /** Inherit types from superclass. */
  typedef typename Superclass::TimeStepType  TimeStepType;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType     FixedImageType;
  typedef typename Superclass::FixedImagePointer  FixedImagePointer;

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType     MovingImageType;
  typedef typename Superclass::MovingImagePointer  MovingImagePointer;
  
   /** Mask image type. */
  typedef TMask MaskImageType;
  typedef typename MaskImageType::Pointer         MaskImagePointer;
  typedef typename MaskImageType::PixelType       MaskPixelType;
  
  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType 
  DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer  
  DeformationFieldPointer;

  /** FiniteDifferenceFunction type. */
  typedef typename Superclass::FiniteDifferenceFunctionType
  FiniteDifferenceFunctionType;

  /** DemonsRegistrationFilterMFunction type. */
  typedef DemonsRegistrationFunctionM<FixedImageType,MovingImageType,
                                     DeformationFieldType,MaskImageType>  DemonsRegistrationFunctionType;
  
  /** Get the metric value. The metric value is the mean square difference 
   * in intensity between the fixed image and transforming moving image 
   * computed over the the overlapping region between the two images. 
   * This is value is only available for the previous iteration and 
   * NOT the current iteration. */
  virtual double GetMetric() const;

  /** Switch between using the fixed image and moving image gradient
   * for computing the deformation field updates. */
  itkSetMacro( UseMovingImageGradient, bool );
  itkGetConstMacro( UseMovingImageGradient, bool );
  itkBooleanMacro( UseMovingImageGradient );

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void SetIntensityDifferenceThreshold(double);
  virtual double GetIntensityDifferenceThreshold() const;
  
protected:
  DemonsRegistrationFilterM();
  ~DemonsRegistrationFilterM() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize the state of filter and equation before each iteration. */
  virtual void InitializeIteration();

  /** Apply update. */
#if ( ITK_VERSION_MAJOR > 3 ) 
  virtual void ApplyUpdate(const TimeStepType &dt);  
#else  
  virtual void ApplyUpdate(TimeStepType dt);
#endif  

private:
  DemonsRegistrationFilterM(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool  m_UseMovingImageGradient;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDemonsRegistrationFilterM.txx"
#endif

#endif
