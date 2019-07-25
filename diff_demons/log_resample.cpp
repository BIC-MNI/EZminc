/* ----------------------------- MNI Header -----------------------------------
@NAME       :  log_resample
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <valarray>
#include <cmath>
#include <limits>
#include <unistd.h>

#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWarpImageFilter.h>

#if ITK_VERSION_MAJOR >= 4
#include <itkExponentialDisplacementFieldImageFilter.h>
#include <itkMultiplyImageFilter.h>
#else
#include <itkMincGeneralTransform.h>
#include <itkExponentialDeformationFieldImageFilter.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <time_stamp.h>    // for creating minc style history entry
#include <itkMincImageIOFactory.h>
#include <itkMincImageIO.h>
#include <itkMincHelpers.h>
#endif

#include <unistd.h>
#include <getopt.h>


#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkImageIOBase.h>

#include "mincUtils.h"

typedef itk::Vector< double, 3 >    VectorType;
typedef itk::Image< VectorType, 3 > VectorImageType;

typedef itk::ImageIOBase          IOBase;
typedef itk::SmartPointer<IOBase> IOBasePointer;

typedef itk::ImageBase<3>           Image3DBase;
typedef itk::Image<float,3>         Float3DImage;
typedef itk::Image<int,3>           Int3DImage;
typedef itk::Image<short,3>         Short3DImage;
typedef itk::Image<unsigned char,3> Byte3DImage;


typedef itk::BSplineInterpolateImageFunction< Float3DImage, double, double >  InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction< Int3DImage, double >    NNInterpolatorType;

typedef itk::ImageFileReader< VectorImageType >  VectorReaderType;

typedef itk::ImageFileReader< Float3DImage > ReaderType;
typedef itk::ImageFileWriter< Float3DImage > WriterType;

#if ITK_VERSION_MAJOR >= 4
typedef itk::ExponentialDisplacementFieldImageFilter< VectorImageType, VectorImageType > ExponentiatorFilterType;
typedef itk::MultiplyImageFilter<VectorImageType, VectorImageType, VectorImageType> MultiplicatorFilterType;
#else    
typedef itk::ExponentialDeformationFieldImageFilter< VectorImageType, VectorImageType > ExponentiatorFilterType;
typedef itk::MultiplyByConstantImageFilter< VectorImageType, double, VectorImageType > MultiplicatorFilterType;
#endif



using namespace  std;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" <input> <output.mnc> " << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--like <example> (default behaviour analogous to use_input_sampling)"<<std::endl
    << "--log_transform <log domain vector field> "<<std::endl
    << "--order <n> spline order, default 2 "<<std::endl
    << "--uniformize <step> - will make a volume with uniform step size and no direction cosines" << std::endl
    << "--invert_transform  - apply inverted transform"<<std::endl
    << "--factor <f>  - multiply by factor"<<std::endl
    << "--labels  - for resampling label data"<<std::endl
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
  typename T::PointType org;
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
  flt->SetOutputSize(size);
  flt->SetOutputOrigin(org);
  flt->SetOutputSpacing(spc);
}


template<class Image,class ImageOut,class Interpolator> 
void resample_image(
   IOBase* base,  
   const std::string& output_f,
   VectorImageType::Pointer input_def_field,
   const std::string& like_f,
   double uniformize,
   const char* history,
   bool store_float,
   bool store_short,
   bool store_byte,
   Interpolator* interpolator)
{

  typedef itk::ImageFileReader< Image > ReaderType;
  typedef itk::ImageFileWriter< ImageOut > WriterType;
  typedef itk::WarpImageFilter < Image, ImageOut, VectorImageType >  WarperType;
  
  typename Image::Pointer in=0;
  {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetImageIO(base);
    reader->SetFileName(base->GetFileName());
    reader->Update();
    in=reader->GetOutput();
    in->DisconnectPipeline();
  }

  //creating the interpolator
  typename WarperType::Pointer filter = WarperType::New();

  filter->SetInterpolator( interpolator );
  filter->SetEdgePaddingValue( 0 );
#if ITK_VERSION_MAJOR >= 4
  filter->SetDisplacementField(input_def_field);
#else
  filter->SetDeformationField(input_def_field);
#endif
  filter->SetNumberOfThreads(1);
  
  if(!like_f.empty())
  {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();
    if(uniformize!=0.0)
    {
      generate_uniform_sampling<WarperType,Image>(filter,reader->GetOutput(),uniformize);
    } else {
      filter->SetOutputParametersFromImage(reader->GetOutput());
      filter->SetOutputDirection(reader->GetOutput()->GetDirection());
    }
  }
  else
  {
    if(uniformize!=0.0)
    {
      generate_uniform_sampling<WarperType,Image>(filter,in,uniformize);
    } else {
      //we are using original sampling
      filter->SetOutputParametersFromImage(in);
      filter->SetOutputDirection(in->GetDirection());
    }
  }
  filter->SetInput(in);
  filter->Update();
  //copy the metadate information, for some reason it is not preserved
  //filter->GetOutput()->SetMetaDataDictionary(reader->GetOutput()->GetMetaDataDictionary());
  typename ImageOut::Pointer out=filter->GetOutput();
  
  //generic file writer
  typename WriterType::Pointer writer = WriterType::New();
  
#ifdef HAVE_MINC4ITK  
  minc::copy_metadata(out,in);
  minc::append_history(out,history);
  
  if(store_float)
  {
    minc::set_minc_storage_type(out,NC_FLOAT,true);
  } else if(store_short) {
    minc::set_minc_storage_type(out,NC_SHORT,true);
  } else if(store_byte) {
    minc::set_minc_storage_type(out,NC_BYTE,false);
  }
#else
  mincify(out,history,store_byte?typeid(unsigned char).name():store_short?typeid(short).name():typeid(float).name(),in);
#endif

  writer->SetFileName(output_f.c_str());
  writer->SetInput( out );
  writer->Update();
}



int main (int argc, char **argv)
{
  int store_float=0;
  int store_short=0;
  int store_byte=0;  

  int verbose=0, clobber=0;
  int order=2;
  std::string like_f,xfm_f,output_f,input_f;
  double uniformize=0.0;
  int invert=0;
  int labels=0;
  double factor=1.0;
  std::string history=minc_timestamp(argc,argv);
  
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"like",    required_argument, 0,       'l'},
    {"log_transform",    required_argument, 0, 't'},
    {"factor",    required_argument, 0, 'f'},
    {"order",    required_argument, 0, 'o'},
    {"uniformize",    required_argument, 0, 'u'},
    {"invert_transform", no_argument, &invert, 1},
    {"labels",  no_argument, &labels, 1},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {0, 0, 0, 0}
    };
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "vqcl:t:o:u:f:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
			{
			case 0:
				break;
			case 'v':
				cout << "Version: 0.1" << endl;
				return 0;
      case 'l':
        like_f=optarg; break;
      case 't':
        xfm_f=optarg; break;
      case 'o':
        order=atoi(optarg);break;
      case 'u':
        uniformize=atof(optarg);break;
      case 'f':
        factor=atof(optarg);break;
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
  
  if(xfm_f.empty())
	{
		std::cerr << "Refusing to work without deformation field!"<<std::endl;
		return 1;
	}
  
	try
  {
#ifdef HAVE_MINC4ITK
    itk::RegisterMincIO();
#endif

    IOBasePointer io = itk::ImageIOFactory::CreateImageIO(input_f.c_str(), itk::ImageIOFactory::ReadMode );
    
    if(!io)
      throw itk::ExceptionObject("Unsupported image file type");
    
    io->SetFileName(input_f.c_str());
    io->ReadImageInformation();

    size_t nd = io->GetNumberOfDimensions();
    size_t nc = io->GetNumberOfComponents();
    
    if(nc!=1) //not supported right now
      throw itk::ExceptionObject("Currently only 1 component images are supported");
    
    if(nd!=3) //not supported right now
      throw itk::ExceptionObject("Currently only 3D images are supported");
    

		VectorImageType::Pointer input_vector_field = 0;
		VectorImageType::Pointer input_def_field = 0;
		
		{
			VectorReaderType::Pointer reader = VectorReaderType::New();
			reader->SetFileName ( xfm_f.c_str() );

			reader->Update();

			input_vector_field = reader->GetOutput();
			input_vector_field->DisconnectPipeline();
		}
		
		{
			ExponentiatorFilterType::Pointer exponentiator = ExponentiatorFilterType::New();
      MultiplicatorFilterType::Pointer multiplicator = MultiplicatorFilterType::New();
      
      if(factor!=1.0)
      {
        multiplicator->SetConstant(factor);
        multiplicator->SetInput(input_vector_field);
        multiplicator->Update();
        exponentiator->SetInput(multiplicator->GetOutput());
      } else 
        exponentiator->SetInput(input_vector_field);
			
			if(!invert) //need to calculate inverse by default!
				exponentiator->ComputeInverseOn();

			exponentiator->Update();
			
			input_def_field = exponentiator->GetOutput();
			input_def_field->DisconnectPipeline();
		}

    if(labels)
    {
      
      NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
      
       if(store_byte)
         resample_image<Int3DImage,Byte3DImage,NNInterpolatorType>(io,output_f,input_def_field,like_f,uniformize,history.c_str(),store_float,store_short,store_byte,interpolator);
       else if(store_short)
         resample_image<Int3DImage,Short3DImage,NNInterpolatorType>(io,output_f,input_def_field,like_f,uniformize,history.c_str(),store_float,store_short,store_byte,interpolator);
       else
         resample_image<Int3DImage,Int3DImage,NNInterpolatorType>(io,output_f,input_def_field,like_f,uniformize,history.c_str(),store_float,store_short,store_byte,interpolator);
      //resample_image<Int3DImage,Int3DImage,NNInterpolatorType>(io,output_f,input_def_field,like_f,uniformize,history,store_float,store_short,store_byte,interpolator);
    } else { 
      InterpolatorType::Pointer interpolator = InterpolatorType::New();
      interpolator->SetSplineOrder(order);
      resample_image<Float3DImage,Float3DImage,InterpolatorType>(io,output_f,input_def_field,like_f,uniformize,history.c_str(),store_float,store_short,store_byte,interpolator);
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
