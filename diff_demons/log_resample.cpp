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
//#include <itkAffineTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <minc_general_transform.h>
#include <itkWarpImageFilter.h>
#include <itkExponentialDeformationFieldImageFilter.h>


#include <unistd.h>
#include <getopt.h>
#include <time_stamp.h>    // for creating minc style history entry


#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>

#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "minc_helpers.h"

//typedef itk::MincImageIO ImageIOType;
typedef itk::Vector< double, 3 >    VectorType;
typedef itk::Image< VectorType, 3 > VectorImageType;
//typedef itk::Image< float, 3 > ImageType;
typedef minc::image3d ImageType;


typedef itk::BSplineInterpolateImageFunction< ImageType, double, double >  InterpolatorType;
//typedef itk::ResampleImageFilter<minc::image3d, minc::image3d> FilterType;
typedef itk::ImageFileReader< VectorImageType >  VectorReaderType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;
typedef itk::ExponentialDeformationFieldImageFilter< VectorImageType, VectorImageType > ExponentiatorFilterType;
//typedef itk::DisplacementFieldCompositionFilter< VectorImageType, VectorImageType >     CompositionFilterType;
typedef itk::WarpImageFilter < ImageType, ImageType, VectorImageType >  WarperType;

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
    << "--invert_transform  - apply inverted transform"<<std::endl;
}

void generate_uniform_sampling(WarperType* flt, const ImageType* img,double step)
{
  //obtain physical coordinats of all image corners
  ImageType::RegionType r=img->GetLargestPossibleRegion();
  std::vector<double> corner[3];
  for(int i=0;i<8;i++)
  {
    ImageType::IndexType idx;
    ImageType::PointType c;
    idx[0]=r.GetIndex()[0]+r.GetSize()[0]*(i%2);
    idx[1]=r.GetIndex()[1]+r.GetSize()[1]*((i/2)%2);
    idx[2]=r.GetIndex()[2]+r.GetSize()[2]*((i/4)%2);
    img->TransformIndexToPhysicalPoint(idx,c);
    for(int j=0;j<3;j++)
      corner[j].push_back(c[j]);
  }
  ImageType::IndexType start;
  WarperType::SizeType size;
  WarperType::PointType org;
  ImageType::SpacingType spc;
  spc.Fill(step);
  for(int j=0;j<3;j++)
  {
    std::sort(corner[j].begin(),corner[j].end());
    size[j]=ceil((corner[j][7]-corner[j][0])/step);
    org[j]=corner[j][0];
  }
  ImageType::DirectionType identity;
  identity.SetIdentity();
  flt->SetOutputDirection(identity);
  start.Fill(0);
  flt->SetOutputStartIndex(start);
  flt->SetOutputSize(size);
  flt->SetOutputOrigin(org);
  flt->SetOutputSpacing(spc);
}

int main (int argc, char **argv)
{
  int verbose=0, clobber=0,skip_grid=0;
  int order=2;
  std::string like_f,xfm_f,output_f,input_f;
  double uniformize=0.0;
  int invert=0;
  char *history = time_stamp(argc, argv); 
  
  static struct option long_options[] = {
		{"verbose", no_argument,       &verbose, 1},
		{"quiet",   no_argument,       &verbose, 0},
		{"clobber", no_argument,       &clobber, 1},
		{"like",    required_argument, 0, 'l'},
		{"log_transform",    required_argument, 0, 't'},
    {"order",    required_argument, 0, 'o'},
    {"uniformize",    required_argument, 0, 'u'},
    {"invert_transform", no_argument, &invert, 1},
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
    itk::RegisterMincIO();
		
		ImageType::Pointer in=0;
		{
			ReaderType::Pointer reader = ReaderType::New();
			
			//initializing the reader
			reader->SetFileName(input_f.c_str());
			reader->Update();
			
			in=reader->GetOutput();
			in->DisconnectPipeline();
		}
		
    WarperType::Pointer filter = WarperType::New();
		

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
			exponentiator->SetInput(input_vector_field);
			
			if(!invert) //need to calculate inverse by default!
				exponentiator->ComputeInverseOn();

			exponentiator->Update();
			
			input_def_field = exponentiator->GetOutput();
			input_def_field->DisconnectPipeline();
		}
		
		
    //creating the interpolator
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		interpolator->SetSplineOrder(order);
		filter->SetInterpolator( interpolator );
		filter->SetEdgePaddingValue( 0 );
    filter->SetDeformationField(input_def_field);
    filter->SetNumberOfThreads(1);
    
    if(!like_f.empty())
    {
      itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
      reader->SetFileName(like_f.c_str());
      reader->Update();
      if(uniformize!=0.0)
      {
        generate_uniform_sampling(filter,reader->GetOutput(),uniformize);
      } else {
        filter->SetOutputParametersFromImage(reader->GetOutput());
        filter->SetOutputDirection(reader->GetOutput()->GetDirection());
      }
    }
    else
    {
      if(uniformize!=0.0)
      {
        generate_uniform_sampling(filter,in,uniformize);
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
    ImageType::Pointer out=filter->GetOutput();
    minc::copy_metadata(out,in);
    minc::append_history(out,history);
    free(history);
    
    //generic file writer
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(output_f.c_str());
    writer->SetInput( out );
    writer->Update();
    
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
