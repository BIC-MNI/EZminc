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
      public  itk::InterpolateImageFunction< TInputImage, TCoordRep >
  {
  public:
    ITK_DISALLOW_COPY_AND_ASSIGN(mincVectorBSplineInterpolate);

    using Self = mincVectorBSplineInterpolate;
    using Superclass = itk::InterpolateImageFunction< TInputImage, TCoordRep >;
    using Pointer = itk::SmartPointer< Self >;
    using ConstPointer = itk::SmartPointer< const Self >;
    
    itkTypeMacro(mincVectorBSplineInterpolate, itk::InterpolateImageFunction);
    
    itkNewMacro(Self);  
    
    /** InputImageType typedef support. */
    using OutputType = typename Superclass::OutputType;
    using InputImageType = typename Superclass::InputImageType;
    using InputPixelType = typename Superclass::InputPixelType;
    using RealType =  TCoordRep; // hack?

    static constexpr unsigned int ImageDimension = Superclass::ImageDimension;
    static constexpr unsigned int Dimension = InputPixelType::Dimension;

    using IndexType = typename Superclass::IndexType;
    using SizeType = typename Superclass::SizeType;
    using ContinuousIndexType = typename Superclass::ContinuousIndexType;
    using InternalComputationType = typename ContinuousIndexType::ValueType;
      
    /** Should check if an index is inside the image buffer, however we
    * require that it answers true to use the extrapolation possibility. */
    // virtual bool IsInsideBuffer( const IndexType & ) const
    // { 
    //   return true;
    // }
    
    /** Evaluate the function at a ContinuousIndex position
    *
    * Returns the linearly interpolated image intensity at a
    * specified point position. If the point does not lie within the
    * image buffer a nearest neighbor interpolation is done. */
    virtual OutputType EvaluateAtContinuousIndex( const ContinuousIndexType & index ) const override;
  
    /** Evaluate the function at an index position
    *
    * Simply returns the image value at the
    * specified index position. If the index does not lie within the
    * image buffer a nearest neighbor interpolation is done. */
    virtual OutputType EvaluateAtIndex( const IndexType & index ) const override;
    
    typedef itk::NthElementImageAdaptor<TInputImage, RealType > ImageAdaptorType;
    typedef typename ImageAdaptorType::Pointer ImageAdaptorPointer;
    
    typedef itk::BSplineInterpolateImageFunction
        <ImageAdaptorType, TCoordRep, RealType >  InterpolatorType;
    typedef typename InterpolatorType::Pointer InterpolatorPointer;
    typedef typename InterpolatorType::CovariantVectorType CovariantVectorType;
    
    virtual void SetInputImage( const InputImageType * ptr );
    
    // CovariantVectorType EvaluateDerivative(unsigned int dim,const PointType &point) const
    // {
    //   return _interpolator[dim]->EvaluateDerivative(point);
    // }
    
/*    CovariantVectorType EvaluateDerivative(unsigned int dim,const PointType &point, unsigned int threadID) const
    {
      return _interpolator[dim]->EvaluateDerivative(point,threadID);
    }*/
    
    // CovariantVectorType EvaluateDerivativeAtContinuousIndex(unsigned int dim, const ContinuousIndexType &x) const
    // {
    //   return _interpolator[dim]->EvaluateDerivativeAtContinuousIndex(x);
    // }
    
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
    
#if !defined(ITKV4_COMPATIBILITY)
    SizeType GetRadius() const override
    {
      return _interpolator[0]->GetRadius();
    }
#endif
    
  protected:
    mincVectorBSplineInterpolate();
    virtual ~mincVectorBSplineInterpolate()  override = default;
  
    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const override ;
  
  private:
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
