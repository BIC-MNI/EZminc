#ifndef __mincVariableVectorBSplineInterpolate_h
#define __mincVariableVectorBSplineInterpolate_h

#include "mincVariableVectorInterpolateImageFunction.h"
#include <itkNthElementImageAdaptor.h>
#include <itkBSplineInterpolateImageFunction.h>

namespace minc
{

  /**
  * \class mincVariableVectorBSplineInterpolate
  * \brief Use B-Spline interpolation to interpolate VectorImage
  *
  *
  * \author Vladimir S. Fonov
  *
  * \warning This function work only for VectorImage 
  *
  * \ingroup ImageFunctions ImageInterpolators
  *
  */
  template <class TInputImage, class TCoordRep = double>
  class mincVariableVectorBSplineInterpolate :
      public  VariableVectorInterpolateImageFunction<TInputImage,TCoordRep>
  {
  public:
    /** Standard class typedefs. */
    typedef mincVariableVectorBSplineInterpolate Self;
    typedef itk::ImageFunction<TInputImage,typename itk::NumericTraits<typename TInputImage::PixelType>::RealType, TCoordRep >     Superclass;
    typedef itk::SmartPointer<Self>                                        Pointer;
    typedef itk::SmartPointer<const Self>                                  ConstPointer;
  
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
  
    /** Run-time type information (and related methods). */
    itkTypeMacro(mincVariableVectorBSplineInterpolate,
                 itk::VectorInterpolateImageFunction);
  
    /** InputImageType typedef support. */
    typedef typename Superclass::InputImageType                           InputImageType;
    typedef typename InputImageType::PixelType                            PixelType;
    typedef typename PixelType::ValueType                                 ValueType;
    typedef typename itk::NumericTraits<ValueType>::RealType              RealType;
    typedef typename Superclass::PointType                                PointType;
  
    /** Dimension underlying input image. */
    itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);
  
    /** Index typedef support. */
    typedef typename Superclass::IndexType                               IndexType;
  
    /** ContinuousIndex typedef support. */
    typedef typename Superclass::ContinuousIndexType                     ContinuousIndexType;
  
    /** Output type is Vector<double,Dimension> */
    typedef typename Superclass::OutputType                              OutputType;
    
    /** Should check if an index is inside the image buffer, however we
    * require that it answers true to use the extrapolation possibility. */
    virtual bool IsInsideBuffer( const IndexType & ) const
    { 
      return true;
    }
  
    /** Should check if a point is inside the image buffer, however we
    * require that it answers true to use the extrapolation possibility. */
    virtual bool IsInsideBuffer( const PointType & ) const
    {
      return true;
    }
  
    /** Should check if a continuous index is inside the image buffer, however we
    * require that it answers true to use the extrapolation possibility. */
    virtual bool IsInsideBuffer( const ContinuousIndexType & ) const
    {
      return true;
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
  
    /** Evaluate the function at a ContinuousIndex position
    *
    * Returns the linearly interpolated image intensity at a
    * specified point position. If the point does not lie within the
    * image buffer a nearest neighbor interpolation is done. */
    virtual OutputType EvaluateAtContinuousIndex(
      const ContinuousIndexType & index ) const;
  
    /** Evaluate the function at an index position
    *
    * Simply returns the image value at the
    * specified index position. If the index does not lie within the
    * image buffer a nearest neighbor interpolation is done. */
    virtual OutputType EvaluateAtIndex( const IndexType & index ) const;
    
    typedef itk::NthElementImageAdaptor<TInputImage, RealType > ImageAdaptorType;
    typedef typename ImageAdaptorType::Pointer ImageAdaptorPointer;
    
    typedef itk::BSplineInterpolateImageFunction
        <ImageAdaptorType, TCoordRep, RealType >  InterpolatorType;
    typedef typename InterpolatorType::Pointer InterpolatorPointer;
    typedef typename InterpolatorType::CovariantVectorType CovariantVectorType;
    
    virtual void SetInputImage( const InputImageType * ptr );
    
    CovariantVectorType EvaluateDerivative(unsigned int dim,const PointType &point) const
    {
      return _interpolator[dim]->EvaluateDerivative(point);
    }
    
/*    CovariantVectorType EvaluateDerivative(unsigned int dim,const PointType &point, unsigned int threadID) const
    {
      return _interpolator[dim]->EvaluateDerivative(point,threadID);
    }*/
    
    CovariantVectorType EvaluateDerivativeAtContinuousIndex(unsigned int dim, const ContinuousIndexType &x) const
    {
      return _interpolator[dim]->EvaluateDerivativeAtContinuousIndex(x);
    }
    
/*    CovariantVectorType EvaluateDerivativeAtContinuousIndex(unsigned int dim, const ContinuousIndexType &x, unsigned int threadID) const
    {
      return _interpolator[dim]->EvaluateDerivativeAtContinuousIndex(x,threadID);
    }*/
    
    void SetSplineOrder (unsigned int SplineOrder)
    {
      /*for(int i=0;i<Dimension;i++)
        _interpolator[i]->SetSplineOrder(SplineOrder);*/
      m_Order=SplineOrder;
    }
    
    unsigned int GetSplineOrder () const 
    {
      return m_Order; //_interpolator[0]->GetSplineOrder();
    }
    
    
  protected:
    mincVariableVectorBSplineInterpolate();
    virtual ~mincVariableVectorBSplineInterpolate() {}
  
    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;
  
  private:
    mincVariableVectorBSplineInterpolate(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
  
    /** Number of neighbors used in the interpolation */
    unsigned long  m_Order;
    
    std::vector<ImageAdaptorPointer> _adaptor;
    std::vector<InterpolatorPointer> _interpolator;
  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "mincVariableVectorBSplineInterpolate.txx"
#endif

#endif
