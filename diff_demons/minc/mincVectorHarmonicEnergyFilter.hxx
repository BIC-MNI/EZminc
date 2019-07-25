#ifndef __itkVectorHarmonicEnergyFilter_hxx
#define __itkVectorHarmonicEnergyFilter_hxx

#include "mincVectorHarmonicEnergyFilter.h"

#include <itkImageRegionIterator.h>
#include <itkProgressReporter.h>
#include <itkMath.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_det.h>

namespace itk
{

template <typename TInputImage, typename TOutputImage>
VectorHarmonicEnergyFilter<TInputImage,  TOutputImage>
::VectorHarmonicEnergyFilter()
{
  //m_RequestedNumberOfThreads = this->GetNumberOfThreads();
  m_RequestedNumberOfThreads=1; //right now BSplineInterpolator is not thread safe :(
  m_Compat=false;
  _interpolate=VectorBSplineInterpolateType::New();
}


template< typename TInputImage,  typename TOutputImage >
void
VectorHarmonicEnergyFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  Superclass::BeforeThreadedGenerateData();
  _interpolate->SetInputImage(this->GetInput());
}

template< typename TInputImage,  typename TOutputImage >
void
VectorHarmonicEnergyFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       itk::ThreadIdType  threadId)
{
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    m_DerivativeWeights[i] = 1.0 / fabs(static_cast<double>(this->GetInput()->GetSpacing()[i]));
  }
  
  itk::ImageRegionIterator<TOutputImage> it(this->GetOutput(), outputRegionForThread);
  
  // Support progress methods/callbacks
  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
 
  while ( ! it.IsAtEnd() )
  {
    vnl_matrix_fixed<double,ImageDimension,VectorDimension> J;
    itk::Point< double, ImageDimension > p;
    this->GetOutput()->TransformIndexToPhysicalPoint(it.GetIndex(),p);

   for(unsigned int i=0;i<VectorDimension;i++)
   {
    typename VectorBSplineInterpolateType::CovariantVectorType _dx=
      _interpolate->EvaluateDerivative(i,p);
      
    for(unsigned int j=0;j<ImageDimension;j++)
      J[j][i]=_dx[j]*m_DerivativeWeights[j];
    //J[i][i] += 1.0; //??
   }
   
   it.Set(J.fro_norm());
    
    ++it;
    progress.CompletedPixel();
  }
}


template <typename TInputImage, typename TOutputImage>
void
VectorHarmonicEnergyFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "m_RequestedNumberOfThreads = " << m_RequestedNumberOfThreads
     << std::endl;
}

} // end namespace itk

#endif
