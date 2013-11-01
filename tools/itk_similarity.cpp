#include <iostream>
#include <itkImageMaskSpatialObject.h>
// cost functions
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkMutualInformationHistogramImageToImageMetric.h>
#include <itkMeanSquaresHistogramImageToImageMetric.h>
#include <itkCorrelationCoefficientHistogramImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkMutualInformationImageToImageMetric.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkKappaStatisticImageToImageMetric.h>

#include <itkResampleImageFilter.h>
//#include <itkIdentityTransform.h>
#include <itkAffineTransform.h>

#include <itkImageFileReader.h>

#if ( ITK_VERSION_MAJOR < 4 )
#include "minc_wrappers.h"
#include <itkMincHelpers.h>
using namespace  minc;
#endif


#include <unistd.h>
#include <getopt.h>

using namespace  std;

//! a helper function for file reading
template <class T> typename T::Pointer load_file(const char *file)
{
  typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
  
  reader->SetFileName(file);
  reader->Update();
  
  return reader->GetOutput();
}


void show_usage (const char *name)
{
  std::cerr 
	  << "Usage: "<<name<<" <src> <target>" << endl
    << "--src_mask <source_mask> " << endl
    << "--target_mask <mask> "<< endl
    << "--mi mutual information"<< endl
    << "--nmi normalized mutual information "<<endl
    << "--mmi Mates Mutual information"<<endl
    << "--msq mean squares "<<endl
    << "--ccoeff Correlation coeffecient"<<endl
    << "--kappa Kappa measure of similarity (binary images)"<<endl;
}

int main (int argc, char **argv)
{
  int verbose=1;
  int mi=0,nmi=0,mmi=0,msq=0,ccoeff=0;
  int kappa=0;
  bool binary=false;
  static struct option long_options[] = { 
		{"src_mask",  required_argument, 0, 's'},
    {"target_mask",  required_argument, 0, 't'},
    {"mi",  no_argument, &mi, 1},
    {"nmi",  no_argument, &nmi, 1},
    {"msq",  no_argument, &msq, 1},
    {"mmi",  no_argument, &mmi, 1},
    {"ccoeff",  no_argument, &ccoeff, 1},
    {"kappa",  no_argument, &kappa, 1},
		{0, 0, 0, 0}
		};
    
  std::string src,target,src_mask,target_mask;
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;
      int c = getopt_long_only(argc, argv, "s:t:m", long_options, &option_index);
      if (c == -1) break;
      
      switch (c)
			{
			case 0:
				break;
      case 's':
        src_mask = optarg;
        break;
      case 't':
        target_mask = optarg;
        break;
			case '?':
				/* getopt_long already printed an error message. */
        //cerr<<"Unknown option:"<<optarg<<endl;
			default:
				show_usage (argv[0]);
				return 1;
			}
  }
	if ((argc - optind) < 2) {
		show_usage(argv[0]);
		return 1;
	}
	src  = argv[optind];
	target = argv[optind+1];
  binary=(kappa==1);
	try
  {
    
    typedef itk::ImageMaskSpatialObject< 3 >   MaskType;
    typedef float FloatVoxel;
    typedef unsigned char BinaryVoxel;
    typedef itk::Image<FloatVoxel, 3>          FloatImage;
    typedef itk::Image<BinaryVoxel, 3>         BinaryImage;
    
    typedef itk::MattesMutualInformationImageToImageMetric<FloatImage,FloatImage> MatesMetric;
    typedef itk::ImageToImageMetric< FloatImage, FloatImage> FloatMetricBase;
    typedef itk::ImageToImageMetric< BinaryImage, BinaryImage>   BinaryMetricBase;
    
    if(binary) // treat input images as binary masks
    {
      BinaryMetricBase::Pointer metric;
      if(kappa) {
        itk::KappaStatisticImageToImageMetric<BinaryImage,BinaryImage>::Pointer  
            kmetric=itk::KappaStatisticImageToImageMetric<BinaryImage,BinaryImage>::New();
        kmetric->SetForegroundValue(1);
        metric=kmetric;
        
      }
      metric->ComputeGradientOff();
      
      BinaryImage::Pointer img_src, img_trg;
      img_src=load_file<BinaryImage>(src.c_str());
      img_trg=load_file<BinaryImage>(target.c_str());
      
      MaskType::Pointer  src_mask_obj( MaskType::New()),
      trg_mask_obj( MaskType::New());
      BinaryImage::Pointer    src_mask_img( BinaryImage::New() ),
      trg_mask_img( BinaryImage::New() );

      typedef itk::NearestNeighborInterpolateImageFunction< BinaryImage,double> NearestNeighborInterpolator;
      typedef itk::AffineTransform<double,3> AffineTransform;
      
      AffineTransform::Pointer transform = AffineTransform::New();
      transform->SetIdentity();
      AffineTransform::ParametersType par=transform->GetParameters();
      NearestNeighborInterpolator::Pointer interpolator=NearestNeighborInterpolator::New();
      
      if(!target_mask.empty())
      {
        trg_mask_img=load_file<BinaryImage>( target_mask.c_str());
        trg_mask_obj->SetImage( trg_mask_img );
        metric->SetFixedImageMask( trg_mask_obj );
      }
      
      if(!src_mask.empty())
      {
        src_mask_img=load_file<BinaryImage>( src_mask.c_str() );
        src_mask_obj->SetImage( src_mask_img );
        metric->SetMovingImageMask( src_mask_obj  );
      }
      
      metric->SetMovingImage(img_src);
      metric->SetFixedImage (img_trg);
      metric->SetFixedImageRegion( img_trg->GetLargestPossibleRegion() );
      interpolator->SetInputImage( img_src );
      
      metric->SetTransform( transform );
      metric->SetInterpolator( interpolator );
      metric->Initialize();
      cout.precision(10);
      cout<<metric->GetValue(par)<<endl;
      
    } else {
      FloatMetricBase::Pointer metric;
      if(ccoeff) {
        metric=itk::CorrelationCoefficientHistogramImageToImageMetric<FloatImage,FloatImage>::New();
      } else if(nmi) {
        metric=itk::NormalizedMutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::New();
      } else if(msq) {
        metric=itk::MeanSquaresImageToImageMetric<FloatImage,FloatImage>::New();
      } else if(mmi) {
        MatesMetric::Pointer m=MatesMetric::New();
        m->SetUseAllPixels(true);
        m->SetNumberOfHistogramBins(50);
        
        metric=m;
      } else /* mi */{
        metric=itk::MutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::New();
      }
      metric->ComputeGradientOff();
      
      FloatImage::Pointer img_src,img_trg;
                       
      img_src=load_file<FloatImage>(src.c_str());
      img_trg=load_file<FloatImage>(target.c_str());
      
      
      MaskType::Pointer  src_mask_obj( MaskType::New()),
                         trg_mask_obj( MaskType::New());
      BinaryImage::Pointer    src_mask_img( BinaryImage::New() ),
                         trg_mask_img( BinaryImage::New() );

      typedef itk::NearestNeighborInterpolateImageFunction< FloatImage,double> NearestNeighborInterpolator;
      typedef itk::AffineTransform<double,3> AffineTransform;
      
      AffineTransform::Pointer transform = AffineTransform::New();
      transform->SetIdentity();
      AffineTransform::ParametersType par=transform->GetParameters();
      NearestNeighborInterpolator::Pointer interpolator=NearestNeighborInterpolator::New();
      
      if(!target_mask.empty())
      {
        trg_mask_img=load_file<BinaryImage>( target_mask.c_str());
        trg_mask_obj->SetImage( trg_mask_img );
        metric->SetFixedImageMask( trg_mask_obj	);
      }
      
      if(!src_mask.empty())
      {
        src_mask_img=load_file<BinaryImage>( src_mask.c_str());
        src_mask_obj->SetImage( src_mask_img );
        metric->SetMovingImageMask( src_mask_obj	);
      }
      
      metric->SetMovingImage(img_src);
      metric->SetFixedImage (img_trg);
      metric->SetFixedImageRegion( img_trg->GetLargestPossibleRegion() );
      interpolator->SetInputImage( img_src );
      
      metric->SetTransform( transform );
      metric->SetInterpolator( interpolator );
      metric->Initialize();
      cout.precision(10);
      cout<<metric->GetValue(par)<<endl;
  //    cerr<<metric<<endl;
    }
      
  } 
#if ( ITK_VERSION_MAJOR < 4 )
  catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  }
#endif  
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
  return 0;
}  

