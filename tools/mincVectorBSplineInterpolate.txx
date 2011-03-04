#ifndef __mincVectorBSplineInterpolate_txx
#define __mincVectorBSplineInterpolate_txx

#include "mincVectorBSplineInterpolate.h"

namespace minc
{

/**
 * Define the number of neighbors
 */
/*template<class TInputImage, class TCoordRep>
const unsigned long mincVectorBSplineInterpolate< TInputImage, TCoordRep>::m_Order = 2;*/


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep> 
mincVectorBSplineInterpolate< TInputImage, TCoordRep >::mincVectorBSplineInterpolate()
  :m_Order(2)
{
  for(unsigned int i=0;i<Dimension;i++)
  {
    _adaptor[i]=ImageAdaptorType::New();
    _interpolator[i]=InterpolatorType::New();
  }
}

template<class TInputImage, class TCoordRep>
void mincVectorBSplineInterpolate< TInputImage, TCoordRep >
::SetInputImage( const InputImageType * ptr )
{
  Superclass::SetInputImage(ptr);
  for(unsigned int i=0;i<Dimension;i++)
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
void mincVectorBSplineInterpolate< TInputImage, TCoordRep >
  ::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename mincVectorBSplineInterpolate< TInputImage, TCoordRep >::OutputType
mincVectorBSplineInterpolate< TInputImage, TCoordRep >
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
  
  for( dim = 0; dim < Dimension; dim++ )
    output[dim]=_interpolator[dim]->EvaluateAtContinuousIndex(insideIndex);

  return ( output );
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename mincVectorBSplineInterpolate< TInputImage, TCoordRep >
::OutputType
mincVectorBSplineInterpolate< TInputImage, TCoordRep >
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

  // Now call the superclass implementation of EvaluateAtIndex
  // since we have ensured that the index lies in the image region
  return this->Superclass::EvaluateAtIndex( insideIndex );
}

} // end namespace itk

#endif
