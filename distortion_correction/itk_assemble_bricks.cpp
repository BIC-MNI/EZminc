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
#include <itkMultiplyImageFilter.h>

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
typedef itk::Image<int,3>           Int3DImage;
typedef itk::Image<short,3>         Short3DImage;
typedef itk::Image<unsigned char,3> Byte3DImage;

typedef itk::BSplineInterpolateImageFunction< Float3DImage, double, double >  InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction< Float3DImage, double >    NNInterpolatorType;
typedef itk::Euler3DTransform<double> RigidTransformType;

using namespace  std;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" <input.tag> <input_brick.mnc> <output.mnc> " << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--like <example> (default behaviour analogous to use_input_sampling)"<<std::endl
    << "--order <n> spline order, default 2 "<<std::endl
    << "--labels - assume that input data is discrete labels, will use Nearest-Neighbour interpolator"<<std::endl
    << "--byte    - store image in byte  voxels minc file"<<std::endl
    << "--short   - store image in short voxels minc file"<<std::endl
    << "--float   - store image in float voxels minc file"<<std::endl
    << "--aa <fwhm> apply anti-aliasing filter to input image before resampling"<<std::endl;
}


template<class Image,class ImageOut,class Interpolator> 
void assemble_bricks(
   const std::string& input_f,
   const std::string& output_f,
   const std::string& like_f,
   const minc::tag_points& tags, 
   const std::vector<int>& tag_labels,
   const std::string& history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator,
   double aa=0.0
                    )
{
  typedef typename itk::ResampleImageFilter<Image, ImageOut> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >   ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >   ImageWriterType;
  typedef typename Image::PixelType InputPixelType;
  typedef typename itk::MaximumImageFilter< ImageOut, ImageOut, ImageOut > MaxFilterType;

  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(input_f.c_str());
  reader->Update();

  typename Image::Pointer input=reader->GetOutput();
  input->DisconnectPipeline();

  //initializing the reader
  typename ResampleFilterType::Pointer filter = ResampleFilterType::New();
  typename MaxFilterType::Pointer  max_filter = MaxFilterType::New();

  //creating coordinate transformation objects
  RigidTransformType::Pointer rigid_transform = RigidTransformType::New();

  //creating the interpolator
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( 0 );

  typename Image::Pointer like=0;
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();

    like=reader->GetOutput();
    like->DisconnectPipeline();
  }

  filter->SetOutputParametersFromImage( like );
  filter->SetOutputDirection( like->GetDirection() );
  filter->SetInput( input );

  typename ImageOut::RegionType region;
  region.SetSize ( like->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( like->GetLargestPossibleRegion().GetIndex() );

  typename ImageOut::Pointer out= ImageOut::New();
  
  out->SetOrigin( like->GetOrigin() );
  out->SetSpacing( like->GetSpacing() );
  out->SetDirection( like->GetDirection() );

  out->SetLargestPossibleRegion( region );
  out->SetBufferedRegion( region );
  out->SetRequestedRegion( region );
  out->Allocate();
  out->FillBuffer(0);

  for(size_t i=0;i<tags.size();i++)
  {
    rigid_transform->SetIdentity();
    
    if(tag_labels[i]>0) 
      rigid_transform->SetRotation(0,0,90.0);

    RigidTransformType::OutputVectorType tr;
    //applying inverse transform
    tr[0]=-tags[i][0];tr[1]=-tags[i][1];tr[2]=-tags[i][2];
    rigid_transform->SetTranslation( tr );

    std::cout<<"\t"<<i<<"\t"<<tr[0]<<","<<tr[1]<<","<<tr[2]<<"\t"<<tag_labels[i]<<std::endl;
    //TODO: force update transform here?
    filter->SetTransform( rigid_transform );

    filter->Update();
    max_filter->SetInput( 0, filter->GetOutput() );
    max_filter->SetInput( 1, out );
    max_filter->Update();

    out=max_filter->GetOutput();
    out->DisconnectPipeline();
  }

  minc::copy_metadata( out , input );
  minc::append_history( out , history );

  //correct dimension order
  if( like.IsNotNull() )
    minc::copy_dimorder( out , like );

  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( output_f.c_str() );

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


template<class Image,class ImageOut,class Interpolator> 
void assemble_bricks_labels (
   const std::string& input_f,
   const std::string& output_f,
   const std::string& like_f,
   const minc::tag_points& tags, 
   const std::vector<int>& tag_labels,
   const std::string& history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator,
   double aa=0.0 )
{
  typedef typename itk::ResampleImageFilter<Float3DImage, Float3DImage> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >                 ImageReaderType;
  typedef typename itk::ImageFileReader<Float3DImage >          FloatImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >              ImageWriterType;
  
  typedef typename itk::MultiplyImageFilter<ImageOut, ImageOut, ImageOut > MultiplyImageFilterType;
  typedef typename itk::MaximumImageFilter< ImageOut, ImageOut, ImageOut > MaxFilterType;

  typedef typename itk::SmoothingRecursiveGaussianImageFilter<Float3DImage,Float3DImage> SmoothingFilterType;
  typedef typename itk::BinaryThresholdImageFilter< Float3DImage, ImageOut > ThresholdFilterType;

  typedef typename Image::PixelType InputPixelType;
  typename FloatImageReaderType::Pointer reader = FloatImageReaderType::New();

  //initializing the reader
  reader->SetFileName(input_f);
  reader->Update();

  typename Float3DImage::Pointer input=reader->GetOutput();
  
  typename ResampleFilterType::Pointer  filter  = ResampleFilterType::New();
  typename MaxFilterType::Pointer  max_filter   = MaxFilterType::New();
  typename MultiplyImageFilterType::Pointer mul_filter = MultiplyImageFilterType::New();
  RigidTransformType::Pointer rigid_transform = RigidTransformType::New();
  
  typename SmoothingFilterType::Pointer smoothing_filter = SmoothingFilterType::New();
  typename ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();

  typename Image::Pointer like=0;

  if(!like_f.empty())
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();

    like=reader->GetOutput();
    like->DisconnectPipeline();
  }

  filter->SetOutputParametersFromImage( like );
  filter->SetOutputDirection( like->GetDirection() );
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( 0.0 );

  if(aa==0.0)
  {
    filter->SetInput( input );
  } else {
    smoothing_filter->SetSigma(aa/(2.0*sqrt(2*log(2.0))));
    smoothing_filter->SetInput( input );
    smoothing_filter->Update();
    
    filter->SetInput( smoothing_filter->GetOutput());
  }
  
  
  //filter->SetNumberOfThreads(1);

  typename ImageOut::RegionType region;
  region.SetSize (like->GetLargestPossibleRegion().GetSize());
  region.SetIndex(like->GetLargestPossibleRegion().GetIndex());

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
    //applying inverse transform
    tr[0]=-tags[i][0];tr[1]=-tags[i][1];tr[2]=-tags[i][2];
    rigid_transform->SetTranslation( tr );

    std::cout<<"\t"<<i<<"\t"<<tr[0]<<","<<tr[1]<<","<<tr[2]<<"\t"<<tag_labels[i]<<std::endl;
    //TODO: force update transform here?
    filter->SetTransform( rigid_transform );
    threshold_filter->SetInput(filter->GetOutput());
    threshold_filter->SetLowerThreshold(0.5);
    threshold_filter->SetUpperThreshold(100);
    threshold_filter->SetInsideValue(1);
    threshold_filter->SetOutsideValue(0);

    mul_filter->SetConstant((i+1));
    mul_filter->SetInput( threshold_filter->GetOutput() );
    max_filter->SetInput( 0, mul_filter->GetOutput() );
    max_filter->SetInput( 1, out );
    max_filter->Update();

    out=max_filter->GetOutput();
    out->DisconnectPipeline();
  }

  //typename ImageOut::Pointer out=filter->GetOutput();
  minc::copy_metadata( out,input );
  minc::append_history( out,history );

  //correct dimension order
  if(like.IsNotNull())
    minc::copy_dimorder(out,like);

  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());

  writer->SetInput( out );
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
  std::string like_f,xfm_f,output_f,brick_f,tag_f;
  
  int invert=0;
  int labels=0;
  std::string history;
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

  if ((argc - optind) < 3) {
    show_usage(argv[0]);
    return 1;
  }

  tag_f   =argv[optind];
  brick_f =argv[optind+1];
  output_f=argv[optind+2];

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
    minc::tag_points tags;
    std::vector<int> tag_labels;

    minc::read_tags(tags,tag_labels,tag_f.c_str());

    if(labels)
    {
      if(order==0 || !order_was_set)
      {
        NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();

        if(store_byte)
          assemble_bricks_labels<Int3DImage,Byte3DImage,NNInterpolatorType> ( brick_f,output_f,like_f,
                                                                      tags,tag_labels,history,
                                                                      store_float,store_short,store_byte,
                                                                      interpolator,
                                                                      aa_fwhm
                                                                            );
        else if(store_short)
          assemble_bricks_labels<Int3DImage,Short3DImage,NNInterpolatorType>( brick_f,output_f,like_f,
                                                                      tags,tag_labels,history,
                                                                      store_float,store_short,store_byte,
                                                                      interpolator,
                                                                      aa_fwhm
                                                                            );
        else
          assemble_bricks_labels<Int3DImage,Int3DImage,NNInterpolatorType>  ( brick_f, output_f, like_f,
                                                                      tags, tag_labels, history,
                                                                      store_float, store_short, store_byte,
                                                                      interpolator,
                                                                      aa_fwhm
                                                                            );
      } else {
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetSplineOrder(order);
        
        if(store_byte)
          assemble_bricks_labels<Int3DImage,Byte3DImage,InterpolatorType> ( brick_f,output_f,like_f,
                                                                      tags,tag_labels,history,
                                                                      store_float,store_short,store_byte,
                                                                      interpolator,
                                                                      aa_fwhm
                                                                            );
        else if(store_short)
          assemble_bricks_labels<Int3DImage,Short3DImage,InterpolatorType>( brick_f,output_f,like_f,
                                                                      tags,tag_labels,history,
                                                                      store_float,store_short,store_byte,
                                                                      interpolator,
                                                                      aa_fwhm
                                                                            );
        else
          assemble_bricks_labels<Int3DImage,Int3DImage,InterpolatorType>  ( brick_f, output_f, like_f,
                                                                      tags, tag_labels, history,
                                                                      store_float, store_short, store_byte,
                                                                      interpolator,
                                                                      aa_fwhm
                                                                            );
      }
    } else {
      InterpolatorType::Pointer interpolator = InterpolatorType::New();
      interpolator->SetSplineOrder(order);
      assemble_bricks<Float3DImage,Float3DImage,InterpolatorType>( brick_f, output_f, like_f,
                                                                   tags, tag_labels, history,
                                                                   store_float, store_short, store_byte,
                                                                   interpolator );
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
