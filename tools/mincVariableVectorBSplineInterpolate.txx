#ifndef __mincVariableVectorBSplineInterpolate_txx
#define __mincVariableVectorBSplineInterpolate_txx

#include "mincVariableVectorBSplineInterpolate.h"

namespace minc
{


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep> 
mincVariableVectorBSplineInterpolate< TInputImage, TCoordRep >::mincVariableVectorBSplineInterpolate()
  :m_Order(2)
{
}

template<class TInputImage, class TCoordRep>
void mincVariableVectorBSplineInterpolate< TInputImage, TCoordRep >
::SetInputImage( const InputImageType * ptr )
{
  Superclass::SetInputImage(ptr);
  
  for(unsigned int i=0;i<m_Dimension;i++)
  {
    _adaptor[i]=ImageAdaptorType::New();
    _interpolator[i]=InterpolatorType::New();
  }
  
  for(unsigned int i=0;i<m_Dimension;i++)
  {
    _adaptor[i]->SetImage((InputImageType*)(ptr));//a hack(?)
    _adaptor[i]->SelectNthElement(i);
    _interpolator[i]->SetSplineOrder(m_Order);
    _interpolator[i]->SetInputImage(_adaptor[i]);
  }
}
/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void mincVariableVectorBSplineInterpolate< TInputImage, TCoordRep >
  ::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename mincVariableVectorBSplineInterpolate< TInputImage, TCoordRep >::OutputType
mincVariableVectorBSplineInterpolate< TInputImage, TCoordRep >
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
  
  for( dim = 0; dim < m_Dimension; dim++ )
    output[dim]=_interpolator[dim]->EvaluateAtContinuousIndex(insideIndex);

  return ( output );
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename mincVariableVectorBSplineInterpolate< TInputImage, TCoordRep >
::OutputType
mincVariableVectorBSplineInterpolate< TInputImage, TCoordRep >
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
  for( unsigned int k = 0; k < m_Dimension; k++ )
  {
    output[k] = static_cast<double>( input[k] );
  }
  return ( output );
}

} // end namespace itk

#endif
