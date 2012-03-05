#ifndef __itkDisplacementToVelocityFieldLogFilter_txx
#define __itkDisplacementToVelocityFieldLogFilter_txx

#include "itkDisplacementToVelocityFieldLogFilter.h"

#include "itkProgressReporter.h"

#include "itkImage.h"
#include "itkImageFileWriter.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
DisplacementToVelocityFieldLogFilter< TInputImage, TOutputImage >
::DisplacementToVelocityFieldLogFilter()
{
  m_NumberOfIterations  = 10;
  m_ElapsedIterations   = 0;
  m_SmoothVelocityField = false;
  m_ExpComp       = ExponentialCompositionFilterType::New();
  m_BCHCalculator = BCHFilterType::New();
  m_LeftSmootherX      = GaussianFilterType::New();
  m_LeftSmootherY      = GaussianFilterType::New();
  m_RightSmootherX     = GaussianFilterType::New();
  m_RightSmootherY     = GaussianFilterType::New();

  m_ExpComp->ComputeInverseOn();
  m_ExpComp->SetNumberOfIntegrationSteps (500);

  m_BCHCalculator->SetNumberOfApproximationTerms(3);

  m_LeftSmootherX->SetDirection (0);
  m_LeftSmootherY->SetDirection (1);
  m_RightSmootherX->SetDirection (0);
  m_RightSmootherY->SetDirection (1);

  m_LeftSmootherX->SetOrder(GaussianFilterType::ZeroOrder);
  m_LeftSmootherY->SetOrder(GaussianFilterType::ZeroOrder);
  m_RightSmootherX->SetOrder(GaussianFilterType::ZeroOrder);
  m_RightSmootherY->SetOrder(GaussianFilterType::ZeroOrder);

  m_LeftSmootherX->SetNormalizeAcrossScale(false);
  m_LeftSmootherY->SetNormalizeAcrossScale(false);
  m_RightSmootherX->SetNormalizeAcrossScale(false);
  m_RightSmootherY->SetNormalizeAcrossScale(false);

  m_LeftSmootherX->SetSigma (2.0);
  m_LeftSmootherY->SetSigma (2.0);
  m_RightSmootherX->SetSigma (2.0);
  m_RightSmootherY->SetSigma (2.0);
}

template< class TInputImage, class TOutputImage >
void
DisplacementToVelocityFieldLogFilter< TInputImage, TOutputImage >
::GenerateData()
{
  typedef ImageRegionIterator< InputImageType > IteratorType;

  // intial value = displacement field
  // v_0 = Phi - Id (Following Bossa's notation)
  typename InputImageType::Pointer current = const_cast< InputImageType * >( this->GetInput() );

  m_ElapsedIterations   = 0;

  ProgressReporter progress(this, 0, m_NumberOfIterations);

  for ( unsigned int i = 0; i < m_NumberOfIterations; i++ )
    {
    // delta_n = exp( -v_n ) o Phi - Id (Following Bossa's notation)
    m_ExpComp->SetVelocityField (current);
    m_ExpComp->SetInput ( this->GetInput() );

    m_ExpComp->Update();

    typename InputImageType::Pointer leftField  = m_ExpComp->GetOutput();
    typename InputImageType::Pointer rightField = current;

    // Smoothing helps stabilizing the computation
    // This was not in Bossa's paper
    if ( m_SmoothVelocityField )
      {
      m_LeftSmootherX->SetInput (leftField);
      m_LeftSmootherY->SetInput ( m_LeftSmootherX->GetOutput() );

      leftField = m_LeftSmootherY->GetOutput();

      m_RightSmootherX->SetInput (rightField);
      m_RightSmootherY->SetInput ( m_RightSmootherX->GetOutput() );

      rightField = m_RightSmootherY->GetOutput();
      }

    // Still following Bossa's notation
    // v_n+1 = v_n + delta_n + 0.5 * [v_n, delta_n] + ...
    m_BCHCalculator->SetInput (1, leftField);
    m_BCHCalculator->SetInput (0, rightField);

    m_BCHCalculator->GetOutput()->SetRequestedRegion ( this->GetOutput()->GetRequestedRegion() );

    m_BCHCalculator->Update();

    current = m_BCHCalculator->GetOutput();
    current->DisconnectPipeline();

    this->GraftOutput (current);

    m_ElapsedIterations++;

    progress.CompletedPixel();     // not really a pixel but an iteration
    this->InvokeEvent( IterationEvent() );
    }

  this->GraftOutput (current);
}

} // end namespace itk

#endif
