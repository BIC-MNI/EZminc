/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkVariableVectorInterpolateImageFunction.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __mincVariableVectorInterpolateImageFunction_h
#define __mincVariableVectorInterpolateImageFunction_h

#include "itkImageFunction.h"
#include "itkFixedArray.h"

namespace itk
{

/** \class VariableVectorInterpolateImageFunction
 * \brief Base class for all VectorImage interpolaters.
 *
 * VariableVectorInterpolateImageFunction is the base for all ImageFunctions that
 * interpolates image with vector pixel types. This function outputs
 * a return value of type Vector<double,Dimension>.
 *
 * This class is templated input image type and the coordinate
 * representation type.
 *
 * \warning This hierarchy of functions work only for VectorImage 
 * 
 * \sa InterpolateImageFunction
 * \ingroup ImageFunctions ImageInterpolators
 */
template <class TInputImage, class TCoordRep = double>
class ITK_EXPORT VariableVectorInterpolateImageFunction :
  public ImageFunction<
    TInputImage,
    typename NumericTraits<typename TInputImage::PixelType>::RealType,
    TCoordRep >
{
public:
  /** Extract the vector dimension from the pixel template parameter. */
  itkStaticConstMacro(Dimension, unsigned int,
                      TInputImage::PixelType::Dimension);
  
  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef VariableVectorInterpolateImageFunction Self;
  typedef ImageFunction<TInputImage,
    typename NumericTraits<typename TInputImage::PixelType>::RealType,
    TCoordRep >                          Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(VariableVectorInterpolateImageFunction, ImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType          InputImageType;
  typedef typename InputImageType::PixelType           PixelType;
  typedef typename PixelType::ValueType                ValueType;
  typedef typename NumericTraits<ValueType>::RealType  RealType;

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType      IndexType;
  typedef typename Superclass::IndexValueType IndexValueType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is RealType of TInputImage::PixelType. */
  typedef typename Superclass::OutputType OutputType;

  /** CoordRep typedef support. */
  typedef TCoordRep CoordRepType;
  
  virtual void SetInputImage( const InputImageType * ptr )
  {
    Superclass::SetInputImage(ptr);
    m_Dimension=ptr->GetNumberOfComponentsPerPixel();
  }


  /** Returns the interpolated image intensity at a 
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType Evaluate( const PointType& point ) const
  {
    ContinuousIndexType index;
    this->GetInputImage()->TransformPhysicalPointToContinuousIndex( point, index );
    return ( this->EvaluateAtContinuousIndex( index ) );
  }

  /** Interpolate the image at a continuous index position
   *
   * Returns the interpolated image intensity at a 
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * Subclasses must override this method.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const = 0;

  /** Interpolate the image at an index position.
   * Simply returns the image value at the
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtIndex( const IndexType & index ) const
  {
    OutputType output;
    PixelType input = this->GetInputImage()->GetPixel( index );
    for( unsigned int k = 0; k < m_Dimension; k++ )
    {
      output[k] = static_cast<double>( input[k] );
    }
    return ( output );
  }

protected:
  VariableVectorInterpolateImageFunction():m_Dimension(0) {}
  
  ~VariableVectorInterpolateImageFunction() {}
  
  void PrintSelf(std::ostream& os, Indent indent) const
    { Superclass::PrintSelf( os, indent ); }

  unsigned int m_Dimension;

private:
  VariableVectorInterpolateImageFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#endif
