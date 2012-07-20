#ifndef __mincVectorBSplineInterpolate_h
#define __mincVectorBSplineInterpolate_h

#include <itkVectorInterpolateImageFunction.h>
#include <itkNthElementImageAdaptor.h>
#include <itkBSplineInterpolateImageFunction.h>

namespace minc
{

  /**
  * \class mincVectorBSplineInterpolate
  * \brief Linearly interpolate or NN extrapolate a vector image at
  * specified positions.
  *
  *
  * \author Vladimir S. Fonov
  *
  * \warning This function work only for Vector images. For
  *
  * \ingroup ImageFunctions ImageInterpolators
  *
  */
  template <class TInputImage, class TCoordRep = double>
  class mincVectorBSplineInterpolate :
      public  itk::VectorInterpolateImageFunction<TInputImage,TCoordRep>
  {
  public:
    /** Standard class typedefs. */
    typedef mincVectorBSplineInterpolate Self;
    typedef itk::VectorInterpolateImageFunction<TInputImage,TCoordRep>     Superclass;
    typedef itk::SmartPointer<Self>                                        Pointer;
    typedef itk::SmartPointer<const Self>                                  ConstPointer;
  
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
  
    /** Run-time type information (and related methods). */
    itkTypeMacro(mincVectorBSplineInterpolate,
                 itk::VectorInterpolateImageFunction);
  
    /** InputImageType typedef support. */
    typedef typename Superclass::InputImageType                           InputImageType;
    typedef typename Superclass::PixelType                                PixelType;
    typedef typename Superclass::ValueType                                ValueType;
    typedef typename Superclass::RealType                                 RealType;
  
    typedef typename Superclass::PointType                                PointType;
  
    /** Grab the vector dimension from the superclass. */
    itkStaticConstMacro(Dimension, unsigned int,
                        Superclass::Dimension);
  
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
    mincVectorBSplineInterpolate();
    virtual ~mincVectorBSplineInterpolate() {}
  
    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;
  
  private:
    mincVectorBSplineInterpolate(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
  
    /** Number of neighbors used in the interpolation */
    unsigned long  m_Order;
    ImageAdaptorPointer _adaptor[Dimension];
    InterpolatorPointer _interpolator[Dimension];
  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "mincVectorBSplineInterpolate.txx"
#endif

#endif
