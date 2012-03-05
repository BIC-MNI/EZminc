#ifndef __itkDisplacementToVelocityFieldLogFilter_h_
#define __itkDisplacementToVelocityFieldLogFilter_h_

#include <itkImageToImageFilter.h>
#include "itkDisplacementFieldCompositionFilter.h"

#if ( ITK_VERSION_MAJOR > 3 ) 
#include <itkExponentialDisplacementFieldImageFilter.h>
#else
#include <itkExponentialDeformationFieldImageFilter.h>
#endif

#include "itkVelocityFieldBCHCompositionFilter.h"
#include "itkVelocityFieldExponentialComposedWithDisplacementFieldFilter.h"
#include <itkRecursiveGaussianImageFilter.h>

namespace itk
{
/** \class DisplacementToVelocityFieldLogFilter
 * \brief Compute the logarithm u of a displacement \phi s.t. u = \exp(\phi)
 *
 * This filter uses the algorithm from
 * M. N. Bossa, S. Olmos Gasso."A new algorithm for the computation of the group
 * logarithm of diffeomorphisms". In Proc. MFCA 2008.
 *
 * \author Pierre Fillard, INRIA Paris
 */

template< class TInputImage, class TOutputImage >
class ITK_EXPORT DisplacementToVelocityFieldLogFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  typedef DisplacementToVelocityFieldLogFilter            Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  itkNewMacro  (Self);
  itkTypeMacro (DisplacementToVelocityFieldLogFilter, ImageToImageFilter);

  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  typedef VelocityFieldBCHCompositionFilter< OutputImageType, OutputImageType > BCHFilterType;
  typedef VelocityFieldExponentialComposedWithDisplacementFieldFilter< TOutputImage, TInputImage, TInputImage >
  ExponentialCompositionFilterType;

  /** Gaussian filtering type. */
  typedef itk::RecursiveGaussianImageFilter< TOutputImage, TOutputImage >
  GaussianFilterType;

  itkSetMacro (NumberOfIterations, unsigned int);
  itkGetMacro (NumberOfIterations, unsigned int);

  itkGetMacro (ElapsedIterations, unsigned int);

  /**
   *  Set/Get the number of integration steps when computing the exponential
   *  of the velocity field (default: 500).
   */
  void         SetNumberOfExponentialIntegrationSteps(unsigned int n)
  {
    m_ExpComp->SetNumberOfIntegrationSteps (n);
  }

  unsigned int GetNumberOfExponentialIntegrationSteps(void) const
  {
    return m_ExpComp->GetNumberOfIntegrationSteps ();
  }

  /**
   * Set/Get the number of approximation terms when computing
   * the BCH (default: 3). */
  void         SetNumberOfBCHApproximationTerms(unsigned int n)
  {
    m_BCHCalculator->SetNumberOfApproximationTerms(n);
  }

  unsigned int GetNumberOfBCHApproximationTerms(void) const
  {
    return m_BCHCalculator->GetNumberOfApproximationTerms();
  }

  /**
   *  Smoothing the velocity field improves stability with large
   *  deformations (but decreases accuracy).
   */
  itkSetMacro     (SmoothVelocityField, bool);
  itkGetMacro     (SmoothVelocityField, bool);
  itkBooleanMacro (SmoothVelocityField);

  void SetSigma(double sigma)
  {
    m_LeftSmootherX->SetSigma (sigma);
    m_LeftSmootherY->SetSigma (sigma);
    m_RightSmootherX->SetSigma (sigma);
    m_RightSmootherY->SetSigma (sigma);
  }

  double GetSigma(void) const
  {
    return m_LeftSmootherX->GetSigma();
  }

protected:
  DisplacementToVelocityFieldLogFilter();
  ~DisplacementToVelocityFieldLogFilter(){}

  void GenerateData(void);

private:
  DisplacementToVelocityFieldLogFilter (const Self &);
  void operator=(const Self &);

  unsigned int m_NumberOfIterations;
  unsigned int m_ElapsedIterations;

  typename ExponentialCompositionFilterType::Pointer m_ExpComp;
  typename BCHFilterType::Pointer m_BCHCalculator;
  typename GaussianFilterType::Pointer m_LeftSmootherX;
  typename GaussianFilterType::Pointer m_LeftSmootherY;
  typename GaussianFilterType::Pointer m_RightSmootherX;
  typename GaussianFilterType::Pointer m_RightSmootherY;

  bool m_SmoothVelocityField;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDisplacementToVelocityFieldLogFilter.txx"
#endif

#endif
