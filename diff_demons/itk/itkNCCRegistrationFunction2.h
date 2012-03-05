#ifndef __itkNCCRegistrationFunction2_h
#define __itkNCCRegistrationFunction2_h

#include "itkPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"

namespace itk {

/**
 * \class NCCRegistrationFunction2
 *
 * This class encapsulate the PDE which drives the demons registration 
 * algorithm. It is used by NCCRegistrationFilter to compute the 
 * output deformation field which will map a moving image onto a
 * a fixed image.
 *
 * Non-integer moving image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetMovingImageInterpolator. Note that the input
 * interpolator must derive from baseclass InterpolateImageFunction.
 *
 * This class is templated over the fixed image type, moving image type,
 * and the deformation field type.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa NCCRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template<class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT NCCRegistrationFunction2 : 
  public PDEDeformableRegistrationFunction< TFixedImage,
    TMovingImage, TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef NCCRegistrationFunction2    Self;
  typedef PDEDeformableRegistrationFunction< TFixedImage,
    TMovingImage, TDeformationField > Superclass;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer<const Self>    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( NCCRegistrationFunction2, 
    PDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType     MovingImageType;
  typedef typename Superclass::MovingImagePointer  MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType     FixedImageType;
  typedef typename Superclass::FixedImagePointer  FixedImagePointer;
  typedef typename FixedImageType::IndexType      IndexType;
  typedef typename FixedImageType::SizeType       SizeType;
  typedef typename FixedImageType::SpacingType    SpacingType;
  
  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType    DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer   
    DeformationFieldTypePointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType           PixelType;
  typedef typename Superclass::RadiusType          RadiusType;
  typedef typename Superclass::NeighborhoodType    NeighborhoodType;
  typedef typename Superclass::FloatOffsetType     FloatOffsetType;
  typedef typename Superclass::TimeStepType        TimeStepType;

  /** Interpolator type. */
  typedef double                                                 CoordRepType;
  typedef InterpolateImageFunction<MovingImageType,CoordRepType> InterpolatorType;

  typedef typename InterpolatorType::Pointer         InterpolatorPointer;
  typedef typename InterpolatorType::PointType       PointType;
  typedef LinearInterpolateImageFunction<MovingImageType,CoordRepType>
                                                     DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double,itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer   GradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void SetMovingImageInterpolator( InterpolatorType * ptr )
    { m_MovingImageInterpolator = ptr; }

  /** Get the moving image interpolator. */
  InterpolatorType * GetMovingImageInterpolator(void)
    { return m_MovingImageInterpolator; }

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void * itkNotUsed(GlobalData)) const
    { return m_TimeStep; }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *GetGlobalDataPointer() const
    {
    GlobalDataStruct *global = new GlobalDataStruct();
    return global;
    }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const
    { delete (GlobalDataStruct *) GlobalData;  }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  virtual PixelType  ComputeUpdate(const NeighborhoodType &neighborhood,
                     void *globalData,
                     const FloatOffsetType &offset = FloatOffsetType(0.0));

  /** Set/Get SubtractMean boolean. If true, the sample mean is subtracted 
   * from the sample values in the cross-correlation formula.
   * Default value is false. */
  void SetSubtractMean( bool e) { m_SubtractMean=e;}
  bool GetSubtractMean( ) const { return m_SubtractMean;}

protected:
  NCCRegistrationFunction2();
  ~NCCRegistrationFunction2() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
    {
    FixedImageNeighborhoodIteratorType   m_FixedImageIterator;
    };

private:
  NCCRegistrationFunction2(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Cache fixed image information. */
  SpacingType                     m_FixedImageSpacing;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer       m_FixedImageGradientCalculator;


  /** Function to interpolate the moving image. */
  InterpolatorPointer             m_MovingImageInterpolator;

  /** The global timestep. */
  TimeStepType                    m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double                          m_DenominatorThreshold;

  mutable double                  m_MetricTotal;

  bool                            m_SubtractMean;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNCCRegistrationFunction2.txx"
#endif

#endif
