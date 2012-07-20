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
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>

#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkMincGeneralTransform.h>
#include <itkVectorImage.h>

#include "mincVariableVectorResampleImageFilter.h"
#include "mincVariableVectorBSplineInterpolate.h"

#include <unistd.h>
#include <getopt.h>
#include <time_stamp.h>    // for creating minc style history entry


#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkImageIOBase.h>

#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "itkMincHelpers.h"

typedef itk::ImageBase<3>   Image3DBase;
typedef itk::Image<float,3> Float3DImage;
typedef itk::Image<int,3>   Int3DImage;
typedef itk::Image<short,3> Short3DImage;
typedef itk::Image<unsigned char,3> Byte3DImage;
//typedef itk::Image<itk::Vector<float,3>, 3>  Vector3DImage;
typedef itk::VectorImage<float, 3>  Vector3DImage;


typedef itk::ImageIOBase          IOBase;
typedef itk::SmartPointer<IOBase> IOBasePointer;

typedef itk::BSplineInterpolateImageFunction< Float3DImage, double, double >  InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction< Int3DImage, double >    NNInterpolatorType;
typedef minc::VariableVectorBSplineInterpolate<Vector3DImage,double>      VectorInterpolatorType;

typedef itk::ResampleImageFilter<Float3DImage, Float3DImage> FloatFilterType;
typedef itk::ResampleImageFilter<Int3DImage  , Int3DImage>   IntFilterType;
typedef itk::VariableVectorResampleImageFilter<Vector3DImage , Vector3DImage>  VectorFilterType;

typedef minc::XfmTransform<double,3,3>  TransformType;

using namespace  std;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" <input> <output.mnc> " << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--like <example> (default behaviour analogous to use_input_sampling)"<<std::endl
    << "--transform <xfm_transform> "<<std::endl
    << "--order <n> spline order, default 2 "<<std::endl
    << "--uniformize <step> - will make a volume with uniform step size and no direction cosines" << std::endl
    << "--invert_transform  - apply inverted transform"<<std::endl
    << "--labels - assume that input data is discrete labels, will use Nearest-Neighbour interpolator"<<std::endl
    << "--byte  - store image in byte  voxels minc file"<<std::endl
    << "--short - store image in short voxels minc file"<<std::endl
    << "--float - store image in float voxels minc file"<<std::endl;
}

template<class T,class I> void generate_uniform_sampling(T* flt, const I* img,double step)
{
  //obtain physical coordinats of all image corners
  typename I::RegionType r=img->GetLargestPossibleRegion();
  std::vector<double> corner[3];
  for(int i=0;i<8;i++)
  {
    typename I::IndexType idx;
    typename I::PointType c;
    idx[0]=r.GetIndex()[0]+r.GetSize()[0]*(i%2);
    idx[1]=r.GetIndex()[1]+r.GetSize()[1]*((i/2)%2);
    idx[2]=r.GetIndex()[2]+r.GetSize()[2]*((i/4)%2);
    img->TransformIndexToPhysicalPoint(idx,c);
    for(int j=0;j<3;j++)
      corner[j].push_back(c[j]);
  }
  typename I::IndexType start;
  typename T::SizeType size;
  typename T::OriginPointType org;
  typename I::SpacingType spc;
  spc.Fill(step);
  for(int j=0;j<3;j++)
  {
    std::sort(corner[j].begin(),corner[j].end());
    size[j]=ceil((corner[j][7]-corner[j][0])/step);
    org[j]=corner[j][0];
  }
  Float3DImage::DirectionType identity;
  identity.SetIdentity();
  flt->SetOutputDirection(identity);
  start.Fill(0);
  flt->SetOutputStartIndex(start);
  flt->SetSize(size);
  flt->SetOutputOrigin(org);
  flt->SetOutputSpacing(spc);
}

template<class Image,class ImageOut,class Interpolator> 
void resample_image(
   IOBase* base,  
   const std::string& output_f,
   const std::string& xfm_f,
   const std::string& like_f,
   bool invert,
   double uniformize,
   const char* history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator)
{
  typedef typename itk::ResampleImageFilter<Image, ImageOut> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >   ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >   ImageWriterType;
  
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  
  //initializing the reader
  reader->SetImageIO(base);
  reader->SetFileName(base->GetFileName());
  //reader->SetFileName(input_f.c_str());
  reader->Update();
  
  typename Image::Pointer in=reader->GetOutput();

  typename ResampleFilterType::Pointer filter  = ResampleFilterType::New();
  
  //creating coordinate transformation objects
  TransformType::Pointer transform = TransformType::New();
  if(!xfm_f.empty())
  {
    //reading a minc style xfm file
    transform->OpenXfm(xfm_f.c_str());
    if(!invert) transform->Invert(); //should be inverted by default to walk through target space
    filter->SetTransform( transform );
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
    if(uniformize!=0.0)
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,reader->GetOutput(),uniformize);
    } else {
      filter->SetOutputParametersFromImage(reader->GetOutput());
      filter->SetOutputDirection(reader->GetOutput()->GetDirection());
    }
    like=reader->GetOutput();
    like->DisconnectPipeline();
  }
  else
  {
    if(uniformize!=0.0)
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,in,uniformize);
    } else {
      //we are using original sampling
      filter->SetOutputParametersFromImage(in);
      filter->SetOutputDirection(in->GetDirection());
    }
  }
  filter->SetInput(in);
  filter->Update();
  typename ImageOut::Pointer out=filter->GetOutput();
  minc::copy_metadata(out,in);
  minc::append_history(out,history);
  
  //correct dimension order
  if(like.IsNotNull())
    minc::copy_dimorder(out,like);
  
  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  
  if(store_float)
  {
    minc::set_minc_storage_type(out,NC_FLOAT,true);
  } else if(store_short) {
    minc::set_minc_storage_type(out,NC_SHORT,true);
  } else if(store_byte) {
    minc::set_minc_storage_type(out,NC_BYTE,false);
  }
  
  writer->SetInput( out );
  writer->Update();
}

template<class Image,class ImageOut,class Interpolator> 
void resample_vector_image(
   IOBase* base,  
   const std::string& output_f,
   const std::string& xfm_f,
   const std::string& like_f,
   bool invert,
   double uniformize,
   const char* history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator)
{
  typedef typename itk::VariableVectorResampleImageFilter<Image, ImageOut> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >   ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >   ImageWriterType;
  
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  
  //initializing the reader
  reader->SetImageIO(base);
  reader->SetFileName(base->GetFileName());
  //reader->SetFileName(input_f.c_str());
  reader->Update();
  
  typename Image::Pointer in=reader->GetOutput();

  typename ResampleFilterType::Pointer filter  = ResampleFilterType::New();
  
  //creating coordinate transformation objects
  TransformType::Pointer transform = TransformType::New();
  if(!xfm_f.empty())
  {
    //reading a minc style xfm file
    transform->OpenXfm(xfm_f.c_str());
    if(!invert) transform->Invert(); //should be inverted by default to walk through target space
    filter->SetTransform( transform );
  }

  //creating the interpolator

  filter->SetInterpolator( interpolator );
  //filter->SetDefaultPixelValue( 0 );
  
  //this is for processing using batch system
  filter->SetNumberOfThreads(1);
  
  typename Image::Pointer like=0;
  if(!like_f.empty())
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();
    if(uniformize!=0.0)
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,reader->GetOutput(),uniformize);
    } else {
      filter->SetOutputSpacing(reader->GetOutput()->GetSpacing());
      filter->SetOutputOrigin(reader->GetOutput()->GetOrigin());
      filter->SetOutputDirection(reader->GetOutput()->GetDirection());
    }
    like=reader->GetOutput();
    like->DisconnectPipeline();
  }
  else
  {
    if(uniformize!=0.0)
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,in,uniformize);
    } else {
      //we are using original sampling
      filter->SetOutputSpacing(in->GetSpacing());
      filter->SetOutputOrigin(in->GetOrigin());
      filter->SetOutputDirection(in->GetDirection());
    }
  }
  filter->SetInput(in);
  filter->Update();
  typename ImageOut::Pointer out=filter->GetOutput();
  minc::copy_metadata(out,in);
  minc::append_history(out,history);
  
  //correct dimension order
  if(like.IsNotNull())
    minc::copy_dimorder(out,like);
  
  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  
  if(store_float)
  {
    minc::set_minc_storage_type(out,NC_FLOAT,true);
  } else if(store_short) {
    minc::set_minc_storage_type(out,NC_SHORT,true);
  } else if(store_byte) {
    minc::set_minc_storage_type(out,NC_BYTE,false);
  }
  
  writer->SetInput( out );
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
  int invert=0;
  int labels=0;
  char *history = time_stamp(argc, argv); 
  
  static struct option long_options[] = {
		{"verbose", no_argument,       &verbose, 1},
		{"quiet",   no_argument,       &verbose, 0},
		{"clobber", no_argument,       &clobber, 1},
		{"like",    required_argument, 0, 'l'},
		{"transform",    required_argument, 0, 't'},
    {"order",    required_argument, 0, 'o'},
    {"uniformize",    required_argument, 0, 'u'},
    {"invert_transform", no_argument, &invert, 1},
    {"labels",no_argument, &labels, 1},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    
    {0, 0, 0, 0}
    };
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "vqcl:t:o:u:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
			{
			case 0:
				break;
			case 'v':
				cout << "Version: 0.8" << endl;
				return 0;
      case 'l':
        like_f=optarg; break;
      case 't':
        xfm_f=optarg; break;
      case 'o':
        order=atoi(optarg);break;
      case 'u':
        uniformize=atof(optarg);break;
			case '?':
				/* getopt_long already printed an error message. */
			default:
				show_usage (argv[0]);
				return 1;
			}
    }

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
    //itk::NthElementImageAdaptor<Vector3DImage,float>::Pointer test=itk::NthElementImageAdaptor<Vector3DImage,float>::New();
    
    if(labels)
      order=0;
    
    itk::RegisterMincIO();
    //try to figure out what we have got
    IOBasePointer io = itk::ImageIOFactory::CreateImageIO(input_f.c_str(), itk::ImageIOFactory::ReadMode );
    
    if(!io)
      throw itk::ExceptionObject("Unsupported image file type");
    
    io->SetFileName(input_f.c_str());
    io->ReadImageInformation();

    size_t nd = io->GetNumberOfDimensions();
    size_t nc = io->GetNumberOfComponents();
    itk::ImageIOBase::IOComponentType  ct = io->GetComponentType();
    itk::ImageIOBase::IOComponentType  oct = ct;
    
//     if(nc!=1) //not supported right now
//       throw itk::ExceptionObject("Currently only 1 component images are supported");
//     
//     if(nd!=3) //not supported right now
//       throw itk::ExceptionObject("Currently only 3D images are supported");

    //std::cout<<"Image IO :"<<io<<std::endl;

    if( nc==1 && nd==3 ) //3D image, simple case
    {
      if(labels)
      {
        //creating the interpolator
        NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
        if(store_byte)
          resample_image<Int3DImage,Byte3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,history,store_float,store_short,store_byte,interpolator);
        else if(store_short)
          resample_image<Int3DImage,Short3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,history,store_float,store_short,store_byte,interpolator);
        else
          resample_image<Int3DImage,Int3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,history,store_float,store_short,store_byte,interpolator);
      } else {
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetSplineOrder(order);
        resample_image<Float3DImage,Float3DImage,InterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,history,store_float,store_short,store_byte,interpolator);
      }
     } else if( nd==3 && nc==3 )  { // deal with this case as vector image 
       VectorInterpolatorType::Pointer interpolator = VectorInterpolatorType::New();
       interpolator->SetSplineOrder(order);
       resample_vector_image<Vector3DImage,Vector3DImage,VectorInterpolatorType>(io,output_f,xfm_f,like_f,invert,uniformize,history,store_float,store_short,store_byte,interpolator);
    } else {
      throw itk::ExceptionObject("This number of dimensions is not supported currently");
    }
    free(history);

    return 0;
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
	return 0;
};
