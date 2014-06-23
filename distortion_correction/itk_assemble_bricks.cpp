/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_resample
@DESCRIPTION:  an example of using spline itnerpolation with MINC xfm transform
@COPYRIGHT  :
              Copyright 2011 Vladimir Fonov, McConnell Brain Imaging Centre, 
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
#include <fstream>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>
#include <algorithm>
#include <stdlib.h>

#include <itkResampleImageFilter.h>
#include <itkEuler3DTransform.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageConstIterator.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>

#include <unistd.h>
#include <getopt.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMaximumImageFilter.h>

#ifdef HAVE_MINC4ITK
#include <time_stamp.h>    // for creating minc style history entry
#include "itkMincHelpers.h"
#include "itkMincGeneralTransform.h"
#include "itkMincImageIO.h"
#include "itkMincHelpers.h"
#else
#include "itk4MincHelpers.h"
#include "itkMINCTransformAdapter.h"
#endif //HAVE_MINC4ITK

typedef itk::ImageBase<3>           Image3DBase;
typedef itk::Image<float,3>         Float3DImage;
typedef itk::Image<short,3>         Int3DImage;


typedef itk::ImageIOBase          IOBase;
typedef itk::SmartPointer<IOBase> IOBasePointer;

typedef itk::BSplineInterpolateImageFunction< Float3DImage, double, double >  InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction< Int3DImage, double >    NNInterpolatorType;
typedef itk::ResampleImageFilter<Float3DImage, Float3DImage> FloatFilterType;
typedef itk::ResampleImageFilter<Int3DImage  , Int3DImage>   IntFilterType;

typedef itk::<double> Euler3DTransform;

using namespace  std;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" <input> <output.mnc> " << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--like <example> (default behaviour analogous to use_input_sampling)"<<std::endl
    << "--order <n> spline order, default 2 "<<std::endl
    << "--labels - assume that input data is discrete labels, will use Nearest-Neighbour interpolator"<<std::endl
    << "--byte  - store image in byte  voxels minc file"<<std::endl
    << "--short - store image in short voxels minc file"<<std::endl
    << "--float - store image in float voxels minc file"<<std::endl
    << "--aa <fwhm> apply anti-aliasing filter to labels before resampling (only usable for order > 0 )"<<std::endl;
}


template<class Image,class ImageOut,class Interpolator> 
void resample_image(
   const std::string& input_f,
   const std::string& output_f,
   const std::string& like_f,
   const tag_points& tags, 
   std::vector<int>& tag_labels,
   const char* history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator )
{
  typedef typename itk::ResampleImageFilter<Image, ImageOut> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >   ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >   ImageWriterType;
  typedef typename Image::PixelType InputPixelType;
  
  typedef typename itk::MaximumImageFilter< ImageOut, ImageOut, ImageOut > MaxFilterType;
  
  typename ImageReaderType::Pointer reader = ImageReaderType::New();

  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(input_f.c_str());
  reader->Update();

  typename Image::Pointer input=reader->GetOutput();
  input->DisconnectPipeline();

  //initializing the reader
  typename ResampleFilterType::Pointer filter  = ResampleFilterType::New();
  MaxFilterType::Pointer max_filter = MaxFilterType::New();

  //creating coordinate transformation objects
  RigidTransformType::Pointer rigid_transform = RigidTransformType::New();

  //creating the interpolator
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( 0 );

  //this is for processing using batch system
  filter->SetNumberOfThreads(1);

  typename Image::Pointer like=0;
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();

    like=reader->GetOutput();
    like->DisconnectPipeline();
  }

  filter->SetOutputParametersFromImage( like->GetOutput() );
  filter->SetOutputDirection( like->GetDirection() );
  filter->SetInput( input );

  typename ImageOut::RegionType region;
  region.SetSize (like->GetSize());
  region.SetIndex(like->GetOutputStartIndex());

  typename ImageOut::Pointer out= ImageOut::New();
  
  out->SetOrigin(like->GetOrigin());
  out->SetSpacing(like->GetSpacing());
  out->SetDirection(like->GetDirection());

  out->SetLargestPossibleRegion(region);
  out->SetBufferedRegion(region);
  out->SetRequestedRegion(region);
  out->Allocate();
  out->FillBuffer(0);

  for(size_t i=0;i<tags.size();i++)
  {
    rigid_transform->SetIdentity();
    if(tag_labels[i]>0) 
      rigid_transform->SetRotation(0,0,90.0);
    RigidTransformType::OutputVectorType tr;
    tr[0]=tags[i][0];tr[1]=tags[i][1];tr[2]=tags[i][2];
    rigid_transform->SetTranslation(tr);

    //TODO: force update transform here?
    filter->SetTransform( rigid_transform );
    
    filter->Update();
    max_filter->SetInput(0,filter->GetOutput);
    max_filter->SetInput(1,out);
    max_filter->Update();
    out=max_filter->Output();
    out->DisconnectPipeline();
  }
  
  minc::copy_metadata(out,in);
  minc::append_history(out,history);

  //correct dimension order
  if(like.IsNotNull())
    minc::copy_dimorder(out,like);
    
  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  
  if( getenv("MINC_COMPRESS") != NULL)
    writer->SetUseCompression( true );
  
#ifdef HAVE_MINC4ITK
  if(store_float)
  {
    minc::set_minc_storage_type(out,NC_FLOAT,true);
  } else if(store_short) {
    minc::set_minc_storage_type(out,NC_SHORT,true);
  } else if(store_byte) {
    minc::set_minc_storage_type(out,NC_BYTE,false);
  }
#else
  if(store_float)
  {
    minc::set_minc_storage_type(out,typeid(float).name());
  } else if(store_short) {
    minc::set_minc_storage_type(out,typeid(unsigned short).name());
  } else if(store_byte) {
    minc::set_minc_storage_type(out,typeid(unsigned char).name());
  }
#endif

  writer->SetInput( out );
  writer->Update();
}

template<class Image,class TmpImage,class ImageOut,class Interpolator> 
void resample_label_image (
   const std::string& input_f,
   const std::string& output_f,
   const std::string& like_f,
   const tag_points& tags, 
   std::vector<int>& labels,
   const char* history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator,
   double aa_fwhm=0.0 )
{
  typedef typename itk::ResampleImageFilter<TmpImage, TmpImage> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >                 ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >              ImageWriterType;
  typedef itk::ImageRegionConstIterator<Image>                  ConstInputImageIteratorType;
  typedef itk::ImageRegionConstIterator<TmpImage>               ConstTmpImageIteratorType;
  typedef itk::ImageRegionIterator<TmpImage>                    TmpImageIteratorType;

  typedef itk::ImageRegionIterator<ImageOut>                    ImageOutIteratorType;
  typedef itk::BinaryThresholdImageFilter<Image,TmpImage>       ThresholdFilterType;
  typedef itk::SmoothingRecursiveGaussianImageFilter<TmpImage,TmpImage>  BlurFilterType;

  typedef typename Image::PixelType InputPixelType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();

  double aa_sigma=aa_fwhm/(2.0*sqrt(2*log(2.0)));
  
  //initializing the reader
  reader->SetImageIO(base);
  reader->SetFileName(base->GetFileName());
  reader->Update();

  typename Image::Pointer in=reader->GetOutput();
  
  if(!label_map.empty())
  {
    for(itk::ImageRegionIterator<Image> it(in, in->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
      std::map<int,int>::const_iterator pos=label_map.find(static_cast<int>(it.Get()));
      if(pos==label_map.end())
        it.Set(0); //set to BG
      else
        it.Set(static_cast<InputPixelType>((*pos).second));
    }
  }

  typename ResampleFilterType::Pointer  filter           = ResampleFilterType::New();
  typename ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
  typename BlurFilterType::Pointer      blur_filter     = BlurFilterType::New();

  //creating coordinate transformation objects
  TransformType::Pointer transform = TransformType::New();
  if(!xfm_f.empty())
  {
#if ( ITK_VERSION_MAJOR < 4 )
    //reading a minc style xfm file
    transform->OpenXfm(xfm_f.c_str());
    if(!invert) transform->Invert(); //should be inverted by default to walk through target space
    filter->SetTransform( transform );
#else
    transform->OpenXfm(xfm_f.c_str());
    if(!invert) transform->Invert(); //should be inverted by default to walk through target space
    filter->SetTransform( transform );
#endif
  }

  //creating the interpolator
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( 0 );
  
  //this is for processing using batch system
  filter->SetNumberOfThreads(1);

  typename Image::Pointer like=0;

  if(!like_f.empty())
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();
    
    like=reader->GetOutput();
    like->DisconnectPipeline();
  }
  
  typename ImageOut::RegionType region;
  region.SetSize (like->GetSize());
  region.SetIndex(like->GetOutputStartIndex());
  
  typename ImageOut::Pointer LabelImage= ImageOut::New();
  LabelImage->SetOrigin(like->GetOrigin());
  LabelImage->SetSpacing(like->GetSpacing());
  LabelImage->SetDirection(like->GetDirection());
  
  LabelImage->SetLargestPossibleRegion(region);
  LabelImage->SetBufferedRegion(region);
  LabelImage->SetRequestedRegion(region);
  LabelImage->Allocate();
  
  
  //typename ImageOut::Pointer out=filter->GetOutput();
  minc::copy_metadata(LabelImage,in);
  minc::append_history(LabelImage,history);
  
  //correct dimension order
  if(like.IsNotNull())
    minc::copy_dimorder(LabelImage,like);
    
  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  
#ifdef HAVE_MINC4ITK
  if(store_float)
  {
    minc::set_minc_storage_type(LabelImage,NC_FLOAT,true);
  } else if(store_short) {
    minc::set_minc_storage_type(LabelImage,NC_SHORT,true);
  } else if(store_byte) {
    minc::set_minc_storage_type(LabelImage,NC_BYTE,false);
  }
#else  
  if(store_float)
  {
    minc::set_minc_storage_type(LabelImage,typeid(float).name());
  } else if(store_short) {
    minc::set_minc_storage_type(LabelImage,typeid(unsigned short).name());
  } else if(store_byte) {
    minc::set_minc_storage_type(LabelImage,typeid(unsigned char).name());
  }
#endif

  writer->SetInput( LabelImage );
  if( getenv("MINC_COMPRESS") != NULL)
    writer->SetUseCompression( true );
  
  writer->Update();
}



int main (int argc, char **argv)
{
  int store_float=0;
  int store_short=0;
  int store_byte=0;
  
  int verbose=0, clobber=0,skip_grid=0;
  int order=2;
  std::string like_f,xfm_f,output_f,input_f;
  double uniformize=0.0;
  double unistep=0.0;
  int normalize=0;
  int invert=0;
  int labels=0;
  std::string history;
  std::string map_f,map_str;
  std::map<int,int> label_map;
  double aa_fwhm=0.0;
  
  
#ifdef HAVE_MINC4ITK
  char *_history = time_stamp(argc, argv); 
  history=_history;
  free(_history);
#else
  history = minc_timestamp(argc,argv);
#endif  
  bool order_was_set=false;
  
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"like",    required_argument, 0, 'l'},
    {"order",   required_argument, 0, 'o'},
    {"labels",  no_argument, &labels, 1},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {"aa",      required_argument,      0,'a'},
    {0, 0, 0, 0}
    };
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "vl:o:a:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
      {
      case 0:
        break;
      case 'v':
        cout << "Version: 0.8" << endl;
        return 0;
      case 'a':
        aa_fwhm=atof(optarg);break;
      case 'l':
        like_f=optarg; break;
      case 'o':
        order=atoi(optarg);
        order_was_set=true;
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
      }
    }

  if(normalize && ! order_was_set) //perform nearest neighbour interpolation in this case
    order=0;

  if ((argc - optind) < 2) {
    show_usage(argv[0]);
    return 1;
  }
  input_f=argv[optind];
  output_f=argv[optind+1];
  
  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }


  try
  {
#if ( ITK_VERSION_MAJOR < 4 )  
    itk::RegisterMincIO();
#endif    
    
    if(labels)
    {
      if(order==0 || !order_was_set)
      {
        //creating the interpolator
        NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
        if(store_byte)
          resample_image<Int3DImage,Byte3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,unistep,normalize,history.c_str(),store_float,store_short,store_byte,interpolator,label_map);
        else if(store_short)
          resample_image<Int3DImage,Short3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,unistep,normalize,history.c_str(),store_float,store_short,store_byte,interpolator,label_map);
        else
          resample_image<Int3DImage,Int3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,unistep,normalize,history.c_str(),store_float,store_short,store_byte,interpolator,label_map);
      } else { //using slow algorithm
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetSplineOrder(order);
        resample_label_image<Int3DImage,Float3DImage,Int3DImage,InterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,unistep,normalize,history.c_str(),store_float,store_short,store_byte,interpolator,label_map,aa_fwhm);
      }
    } else {
      InterpolatorType::Pointer interpolator = InterpolatorType::New();
      interpolator->SetSplineOrder(order);
      resample_image<Float3DImage,Float3DImage,InterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,unistep,normalize,history.c_str(),store_float,store_short,store_byte,interpolator);
    }
    return 0;
  } 
#ifdef HAVE_MINC4ITK
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
};
