/* ----------------------------- MNI Header -----------------------------------
@NAME       : itk_vesselness
@DESCRIPTION: apply vesselness filter to the volume
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
#ifdef HAVE_MINC4ITK
  #include <time_stamp.h>    // for creating minc style history entry
  #include "itkMincImageIOFactory.h"
  #include "itkMincImageIO.h"
  #include "itkMincHelpers.h"
#else
  #include "itk4MincHelpers.h"
#endif

#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include <iostream>
#include <getopt.h>
#include <unistd.h>
#include <itkHessian3DToVesselnessMeasureImageFilter.h>
#include <itkMultiScaleHessianBasedMeasureImageFilter.h>
#include <itkHessianToObjectnessMeasureImageFilter.h>

using namespace  std;
using namespace  minc;

void show_usage (const char *name)
{
  std::cerr 
    << "Usage: " << name << " <input> <output> " << endl
    << "--clobber clobber output files" << endl
    << "--version" << endl
    << "--sigma <d>  default 1 " << endl
    << "--alpha <d> " << endl
    << "--beta <d> " << endl
    << "--gamma <d> " << endl
    << "--bright" << endl
    << "--frangi - use Frangi's filter (2nd order derivative)" << endl
    << "--short  - store image in short file format" << endl
    << "--scales <n> - use multiscale algo" << endl
    << "--sigma2 <f> - second sigma for multiscale" << endl;
}


int main (int argc, char **argv)
{
#if ( ITK_VERSION_MAJOR < 4 )
  char *_history = time_stamp(argc, argv); 
  std::string history=_history;
  free(_history);
#else
  std::string history= minc_timestamp(argc,argv);  
#endif  
  
  int clobber=0;
  int c;
  double sigma=0.0;
  double sigma2=0.0;
  double alpha1=0.5;
  double alpha2=0.5;
  double gamma=5.0;
  double beta=0.5;
  int bright_object=0;
  int frangi=0;
  int scales=0;
  int store_byte=0;
  int store_short=0;
  int store_float=0;
  
  static struct option long_options[] = { 
    {"clobber", no_argument, &clobber, 1},
    {"frangi",  no_argument, &frangi, 1},
    {"short",   no_argument, &store_short, 1},
    {"float",   no_argument, &store_float, 1},
    {"byte",    no_argument, &store_byte, 1},
    {"version", no_argument, 0, 'v'},
    {"sigma",  required_argument, 0, 's'},
    {"sigma2", required_argument, 0, 'S'},
    {"alpha1", required_argument, 0, 'a'},
    {"alpha",  required_argument, 0, 'a'},
    {"alpha2", required_argument, 0, 'b'},
    {"beta",   required_argument, 0, 'b'},
    {"gamma",  required_argument, 0, 'g'},
    {"scales", required_argument, 0, 'm'},
    {"bright", no_argument, &bright_object, 1},
    {0, 0, 0, 0}};

	for (;;)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "cvs:a:b:g:s:S:m:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
		case 's':
			sigma = atof(optarg);
			break;
    case 'a':
      alpha1 = atof(optarg);
      break;
    case 'b':
      beta=alpha2 = atof(optarg);
      break;
    case 'g':
      gamma=atof(optarg);
      break;
    case 'S':
      sigma2=atof(optarg);
      break;
    case 'm':
      scales=atoi(optarg);
      break;
    case 'v':
      cout << "Version: 1.0" << endl;
      return 0;
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<2)
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string input(argv[optind]),output(argv[optind+1]);
  
  if (!clobber && !access (output.c_str (), F_OK))
  {
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
	try
  {
#ifdef HAVE_MINC4ITK
    itk::RegisterMincIO();
#endif
    //double sigma=fwhm/(2.0*sqrt(2*log(2.0)));

    itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
    reader->SetFileName(input.c_str());
    
    typedef itk::HessianRecursiveGaussianImageFilter<  image3d > HessianFilterType;    
    
    typedef itk::Hessian3DToVesselnessMeasureImageFilter< float > VesselnessMeasureFilterType;
              
    typedef itk::HessianToObjectnessMeasureImageFilter< HessianFilterType::OutputImageType,image3d >
          FrangiFilterType;
          
//     typedef itk::HessianToObjectnessMeasureImageFilter<float,3> 
//           ObjectnessFilterType;
    
    typedef itk::MultiScaleHessianBasedMeasureImageFilter< image3d, VesselnessMeasureFilterType, image3d > 
          MultiScaleFilterType;
    
    itk::Image<float, 3>::Pointer out;
    
    if(scales>1)
    {
        MultiScaleFilterType::Pointer filter = MultiScaleFilterType::New();
        
        filter->SetInput(imageReader->GetOutput());
        filter->SetSigmaMin(sigma1);
        filter->SetSigmaMax(sigma2);
        filter->SetNumberOfSigmaSteps(scales);
        
        VesselnessMeasureFilterType* vesselness =
            filter->GetHessianToMeasureFilter();
            
        vesselness->SetScaleObjectnessMeasure(false);
        vesselness->SetBrightObject(bright_object);
        vesselness->SetAlpha1(alpha1);
        vesselness->SetAlpha2(alpha2);
//         vesselness->SetBeta(beta);
//         vesselness->SetGamma(gamma);
//         vesselness->SetObjectDimension(1);        
        out=filter->GetOutput();

    }
    else if(frangi)
    {
      HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
      FrangiFilterType::Pointer frangiFilter = FrangiFilterType::New();

      if(sigma>0.0)
        hessianFilter->SetSigma( sigma );
        
      hessianFilter->SetInput( reader->GetOutput() );
      
      frangiFilter->SetObjectDimension(2);
      frangiFilter->SetBrightObject(bright_object);
      
      frangiFilter->SetInput( hessianFilter->GetOutput() );

      if( alpha1>0.0 )
      {
        frangiFilter->SetAlpha( alpha1 );
      }

      if( alpha2>0.0 )
      {
        frangiFilter->SetBeta( beta );
      }

      if( gamma >0.0)
      {
        frangiFilter->SetGamma( gamma );
      }
      frangiFilter->Update();
      
      out=frangiFilter->GetOutput();
      
    } else {
      HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
      VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();
      hessianFilter->SetInput( reader->GetOutput() );
      
      if( alpha1>0.0 )
      {
        vesselnessFilter->SetAlpha1( alpha1 );
      }

      if( alpha2>0.0 )
      {
        vesselnessFilter->SetAlpha2( alpha2 );
      }

      vesselnessFilter->SetInput( hessianFilter->GetOutput() );
      vesselnessFilter->Update();
      
      out=vesselnessFilter->GetOutput();
    }
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
    minc::append_history(out,history);

    itk::ImageFileWriter< itk::Image<float, 3> >::Pointer writer = itk::ImageFileWriter<itk::Image<float, 3> >::New();
    
    writer->SetFileName(output.c_str());
    writer->SetInput( out );
    writer->Update();
    
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
}
