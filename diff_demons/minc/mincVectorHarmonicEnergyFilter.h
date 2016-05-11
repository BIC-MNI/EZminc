#ifndef __itkVectorHarmonicEnergyFilter_h
#define __itkVectorHarmonicEnergyFilter_h

#include <itkImageToImageFilter.h>
#include "mincVectorBSplineInterpolate.h"
#include <itkVector.h>

namespace itk
{
  template < typename TInputImage,
             typename TOutputImage  >
  class VectorHarmonicEnergyFilter :
      public itk::ImageToImageFilter< TInputImage, TOutputImage >
  {
  public:
    /** Standard class typedefs. */
    typedef VectorHarmonicEnergyFilter      Self;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self>                              Pointer;
    typedef itk::SmartPointer<const Self>                        ConstPointer;
  
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
  
    /** Run-time type information (and related methods) */
    itkTypeMacro(VectorHarmonicEnergyFilter, ImageToImageFilter);
  
    /** Extract some information from the image types.  Dimensionality
    * of the two images is assumed to be the same. */
    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef typename TInputImage::PixelType  InputPixelType;
  
    /** Image typedef support */
    typedef TInputImage                       InputImageType;
    typedef TOutputImage                      OutputImageType;
    typedef typename InputImageType::Pointer  InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;
  
    /** The dimensionality of the input and output images. */
    itkStaticConstMacro(ImageDimension, unsigned int,
                        TOutputImage::ImageDimension);
  
    /** Length of the vector pixel type of the input image. */
    itkStaticConstMacro(VectorDimension, unsigned int,
                        InputPixelType::Dimension);
  
    /** Define the data type and the vector of data type used in calculations. */
    typedef OutputPixelType RealType;
  
    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    
    typedef mincVectorBSplineInterpolate<InputImageType>   VectorBSplineInterpolateType;
    typedef typename VectorBSplineInterpolateType::Pointer VectorBSplineInterpolatePointer;
  
    /** mincblob compatibility mode (subtract 1 from determinant) */
   
    itkSetMacro( Compat, bool );
    itkGetMacro( Compat, bool );
    itkBooleanMacro(Compat);
  
  protected:
    VectorHarmonicEnergyFilter();
    virtual ~VectorHarmonicEnergyFilter() {}
  
    /** Do any necessary casting/copying of the input data.  Input pixel types
      whose value types are not real number types must be cast to real number
      types. */
    void BeforeThreadedGenerateData ();
  
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                              itk::ThreadIdType  threadId );
  
    void PrintSelf(std::ostream& os, itk::Indent indent) const;
  
    typedef typename InputImageType::Superclass ImageBaseType;
  
    /** Get access to the input image casted as real pixel values */
    itkGetConstObjectMacro( RealValuedInputImage, ImageBaseType );
    

 
  private:
    int  m_RequestedNumberOfThreads;
    bool  m_Compat;
    double m_DerivativeWeights[ImageDimension];
    
    VectorBSplineInterpolatePointer _interpolate;
  
    typename ImageBaseType::ConstPointer m_RealValuedInputImage;
  
    VectorHarmonicEnergyFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "mincVectorHarmonicEnergyFilter.hxx"
#endif

#endif
