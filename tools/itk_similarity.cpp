/* ----------------------------- MNI Header -----------------------------------
@NAME       : itk_similarity
@DESCRIPTION: calculate various ITK image similarity metrics
@COPYRIGHT  :
              Copyright 2006 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include <iostream>
#include <itkImageMaskSpatialObject.h>
// cost functions

#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkMutualInformationHistogramImageToImageMetric.h>
#include <itkMutualInformationImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkCorrelationCoefficientHistogramImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkKullbackLeiblerCompareHistogramImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkKappaStatisticImageToImageMetric.h>

#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>

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
  
  typename T::Pointer r=reader->GetOutput();
  return r;
}


void show_usage (const char *name)
{
  std::cerr 
    << "Usage: "<<name<<" <src> <target>" << endl
    << "--verbose print metric name"      << endl
    << "--src_mask <source_mask> "        << endl
    << "--target_mask <mask> "            << endl
    << "--mi mutual information"<< endl
    << "--nmi normalized mutual information "<<endl
    << "--mmi Mates Mutual information"<<endl
    //<< "--kl Kullback Leibler histogram similarity"<<endl
    << "--msq mean squares "<<endl
    << "--cc correlation coeffecient"<<endl
    << "--ncc Normalized correlation coeffecient"<<endl
    << "--kappa Kappa measure of similarity (binary images)"<<endl
    << "--bins <n> set number of histogram bins (per dimension)"<<endl;
}

int main (int argc, char **argv)
{
  int verbose=0;
  int mi=0,nmi=0,mmi=0,msq=0,ccoeff=0,kl=0,ncc=0;
  
  int  kappa=0;
  bool binary=false;
  int  hist_bins=64;
  
  static struct option long_options[] = { 
    {"src_mask",  required_argument, 0, 's'},
    {"target_mask",  required_argument, 0, 't'},
    {"bins",  required_argument, 0, 'b'},
    {"mi",   no_argument, &mi, 1},
    {"nmi",  no_argument, &nmi, 1},
    {"msq",  no_argument, &msq, 1},
    {"mmi",  no_argument, &mmi, 1},
    {"ccoeff",  no_argument, &ccoeff, 1},
    {"cc",  no_argument, &ccoeff, 1},
    {"ncc",  no_argument, &ncc, 1},
    {"kl",  no_argument, &kl, 1},
    {"kappa",  no_argument, &kappa, 1},
    {"verbose",  no_argument, &verbose, 1},
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
      case 'b':
        hist_bins=atoi(optarg);
        break;
      case '?':
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
#ifdef HAVE_MINC4ITK
    itk::RegisterMincIO();
#endif

    typedef itk::ImageMaskSpatialObject< 3 >   MaskType;
    typedef float FloatVoxel;
    typedef unsigned char BinaryVoxel;
    typedef itk::Image<FloatVoxel, 3>          FloatImage;
    typedef itk::Image<BinaryVoxel, 3>         BinaryImage;
    
    typedef itk::ImageToImageMetric< FloatImage,   FloatImage>   FloatMetricBase;
    typedef itk::ImageToImageMetric< BinaryImage, BinaryImage>   BinaryMetricBase;

    // sign for metric - to make sure all metrics show higher value when similarity is maximized
    int sign=1;
    
    if(binary) // treat input images as binary masks
    {
      BinaryMetricBase::Pointer metric;
      //if(kappa) {
        itk::KappaStatisticImageToImageMetric<BinaryImage,BinaryImage>::Pointer  
            kmetric=itk::KappaStatisticImageToImageMetric<BinaryImage,BinaryImage>::New();
        kmetric->SetForegroundValue(1);
        metric=kmetric;
        if(verbose) cout<<"kappa=";
      //}
      metric->ComputeGradientOff();
      sign=1;
      
      BinaryImage::Pointer img_src, img_trg;
      img_src=load_file<BinaryImage>(src.c_str());
      img_trg=load_file<BinaryImage>(target.c_str());
      
      MaskType::Pointer  src_mask_obj( MaskType::New()),
      trg_mask_obj( MaskType::New());
      BinaryImage::Pointer    src_mask_img( BinaryImage::New() ),
      trg_mask_img( BinaryImage::New() );

      typedef itk::NearestNeighborInterpolateImageFunction< BinaryImage,double> InterpolatorType;
      typedef itk::IdentityTransform<double,3> TransformType;
      
      TransformType::Pointer transform = TransformType::New();
      TransformType::ParametersType par=transform->GetParameters();
      
      InterpolatorType::Pointer interpolator=InterpolatorType::New();
      
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
      metric->ComputeGradientOff();
      metric->SetUseAllPixels(true);
      
      metric->SetMovingImage(img_src);
      metric->SetFixedImage (img_trg);
      
      metric->SetFixedImageRegion( img_trg->GetLargestPossibleRegion() );
      interpolator->SetInputImage( img_src );
      
      metric->SetTransform( transform );
      metric->SetInterpolator( interpolator );
      metric->Initialize();
      
      cout.precision(10);
      cout<<metric->GetValue(par)*sign<<endl;
      
    } else {
      FloatMetricBase::Pointer metric;
      
      if(ccoeff) {
        itk::CorrelationCoefficientHistogramImageToImageMetric<FloatImage,FloatImage>::Pointer m=
        itk::CorrelationCoefficientHistogramImageToImageMetric<FloatImage,FloatImage>::New();
        sign=1;

        typedef itk::CorrelationCoefficientHistogramImageToImageMetric<FloatImage,FloatImage>::HistogramSizeType HistogramSizeType;
        HistogramSizeType histogramSize;
        histogramSize.SetSize(2);
        histogramSize[0] = hist_bins;
        histogramSize[1] = hist_bins;
        m->SetHistogramSize( histogramSize );

        metric=m;
        if(verbose) cout<<"cc=";
      } else if(nmi) {
        itk::NormalizedMutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::Pointer m=itk::NormalizedMutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::New();
        sign=1;
        
        typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::HistogramSizeType HistogramSizeType;
        HistogramSizeType histogramSize;
        histogramSize.SetSize(2);
        histogramSize[0] = hist_bins;
        histogramSize[1] = hist_bins;
        m->SetHistogramSize( histogramSize );

        metric=m;
        if(verbose) cout<<"nmi=";
      } else if(msq) {
        metric=itk::MeanSquaresImageToImageMetric<FloatImage,FloatImage>::New();
        
      } else if(mmi) {
        itk::MattesMutualInformationImageToImageMetric<FloatImage,FloatImage>::Pointer m=itk::MattesMutualInformationImageToImageMetric<FloatImage,FloatImage>::New();
        sign=-1;
        
        m->SetUseAllPixels(true);
        m->SetNumberOfHistogramBins(hist_bins);
        
        metric=m;
        if(verbose) cout<<"mmi=";
      } else if(kl) {
        itk::KullbackLeiblerCompareHistogramImageToImageMetric<FloatImage,FloatImage>::Pointer m=itk::KullbackLeiblerCompareHistogramImageToImageMetric<FloatImage,FloatImage>::New();

        typedef itk::KullbackLeiblerCompareHistogramImageToImageMetric<FloatImage,FloatImage>::HistogramSizeType HistogramSizeType;
        HistogramSizeType histogramSize;
        histogramSize.SetSize(2);
        histogramSize[0] = hist_bins;
        histogramSize[1] = hist_bins;
        m->SetHistogramSize( histogramSize );

        metric=m;
        if(verbose) cout<<"kl=";
      } else if(ncc) {
        itk::NormalizedCorrelationImageToImageMetric<FloatImage,FloatImage>::Pointer m=itk::NormalizedCorrelationImageToImageMetric<FloatImage,FloatImage>::New();

        metric=m;
        if(verbose) cout<<"ncc=";
      } else /* mi */{
        itk::MutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::Pointer m=itk::MutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::New();

        typedef itk::MutualInformationHistogramImageToImageMetric<FloatImage,FloatImage>::HistogramSizeType HistogramSizeType;
        HistogramSizeType histogramSize;
        histogramSize.SetSize(2);
        histogramSize[0] = hist_bins;
        histogramSize[1] = hist_bins;
        
        m->SetHistogramSize( histogramSize );
        
        metric=m;
        if(verbose) cout<<"mi=";
      }
      
      metric->ComputeGradientOff();
      metric->SetUseAllPixels(true);
      
      FloatImage::Pointer img_src,img_trg;
                       
      img_src=load_file<FloatImage>(src.c_str());
      img_trg=load_file<FloatImage>(target.c_str());
      
      
      MaskType::Pointer  src_mask_obj( MaskType::New()),
                         trg_mask_obj( MaskType::New());

      BinaryImage::Pointer    src_mask_img( BinaryImage::New() ),
                              trg_mask_img( BinaryImage::New() );

      typedef itk::NearestNeighborInterpolateImageFunction< FloatImage,double> InterpolatorType;
      typedef itk::AffineTransform<double,3> AffineTransform;
      
      typedef itk::IdentityTransform<double,3> TransformType;
      
      TransformType::Pointer transform = TransformType::New();
      transform->SetIdentity();
      TransformType::ParametersType par=transform->GetParameters();
      InterpolatorType::Pointer interpolator=InterpolatorType::New();
      
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
      cout<<metric->GetValue(par)*sign<<endl;
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

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;tab-width 2

