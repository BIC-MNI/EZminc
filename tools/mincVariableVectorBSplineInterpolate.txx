#ifndef __mincVariableVectorBSplineInterpolate_txx
#define __mincVariableVectorBSplineInterpolate_txx

#include "mincVariableVectorBSplineInterpolate.h"

namespace minc
{


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep> 
VariableVectorBSplineInterpolate< TInputImage, TCoordRep >::VariableVectorBSplineInterpolate()
  :m_Order(2)
{
}

template<class TInputImage, class TCoordRep>
void VariableVectorBSplineInterpolate< TInputImage, TCoordRep >
::SetInputImage( const InputImageType * ptr )
{
  Superclass::SetInputImage(ptr);
  
  for(unsigned int i=0;i<this->GetNumberOfComponentsPerPixel();i++)
  {
    _adaptor[i]=ImageAdaptorType::New();
    _caster[i]=ImageCastFilterType::New();
    _interpolator[i]=InterpolatorType::New();
  }
  
  for(unsigned int i=0;i<this->GetNumberOfComponentsPerPixel();i++)
  {
    _adaptor[i]->SetImage((InputImageType*)(ptr));//a hack(?)
    _adaptor[i]->SetExtractComponentIndex(i);
    _caster[i]->SetInputImage(_adaptor[i]);
    //_caster[i]->Update();
    _interpolator[i]->SetSplineOrder(m_Order);
    _interpolator[i]->SetInputImage(_caster[i]->GetOutput());
  }
}
/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void VariableVectorBSplineInterpolate< TInputImage, TCoordRep >
  ::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename VariableVectorBSplineInterpolate< TInputImage, TCoordRep >::OutputType
VariableVectorBSplineInterpolate< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex( const ContinuousIndexType& index) const
{
  unsigned int dim;  // index over dimension

  OutputType output;
  output.Fill( 0.0 );
  ContinuousIndexType insideIndex=index;
  
  for( dim = 0; dim < ImageDimension; dim++ )
  {
    if( insideIndex[dim] <  this->m_StartIndex[dim] )
    {
      insideIndex[dim]=this->m_StartIndex[dim];
    }
    if(insideIndex[dim] >=  this->m_EndIndex[dim])
    {
      insideIndex[dim]=this->m_EndIndex[dim]-1;
    }
  }
  
  for( dim = 0; dim < this->GetNumberOfComponentsPerPixel(); dim++ )
    output[dim]=_interpolator[dim]->EvaluateAtContinuousIndex(insideIndex);

  return ( output );
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename VariableVectorBSplineInterpolate< TInputImage, TCoordRep >
::OutputType
VariableVectorBSplineInterpolate< TInputImage, TCoordRep >
::EvaluateAtIndex( const IndexType & index) const
{
  unsigned int dim;
  IndexType insideIndex=index;
  
  for( dim = 0; dim < ImageDimension; dim++ )
  {
    if( insideIndex[dim] <  this->m_StartIndex[dim] )
    {
      insideIndex[dim]=this->m_StartIndex[dim];
    }
    if(insideIndex[dim] >=  this->m_EndIndex[dim])
    {
      insideIndex[dim]=this->m_EndIndex[dim]-1;
    }
  }

  OutputType output;
  PixelType input = this->GetInputImage()->GetPixel( index );
  for( unsigned int k = 0; k < this->GetNumberOfComponentsPerPixel(); k++ )
  {
    output[k] = static_cast<double>( input[k] );
  }
  return ( output );
}

} // end namespace minc

#endif
