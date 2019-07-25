/* ----------------------------- MNI Header -----------------------------------
@NAME       : itk_diffusion
@DESCRIPTION: apply anisotropi diffusion
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
#include <unistd.h>

#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>

#include <iostream>
#include <getopt.h>

#ifdef HAVE_MINC4ITK
#include <time_stamp.h>    // for creating minc style history entry
#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "itkMincHelpers.h"
#else
#include "itk4MincHelpers.h"
#endif 


typedef itk::GradientAnisotropicDiffusionImageFilter< minc::image3d, minc::image3d>  GradientAnisotropicFilter;
typedef itk::CurvatureAnisotropicDiffusionImageFilter< minc::image3d, minc::image3d>  CurvatureAnisotropicFilter;
typedef itk::AnisotropicDiffusionImageFilter< minc::image3d, minc::image3d>  AnisotropicFilter;

using namespace  std;
using namespace  minc;

void show_usage (const char *name)
{
  std::cerr 
    << "Minc wrapper around itk::CurvatureAnisotropicDiffusionImageFilter and itk::GradientAnisotropicDiffusionImageFilter"<<std::endl
    << "Usage: " << name << " <input> <output> " << endl
    << "--clobber clobber output files" << endl
    << "--curvature for curvature anisotropic diffusion filter"<<endl
    << "--gradient  for gradient anisotropic diffusion filter (default)"<<endl
    << "--step <f> - time step, default 0.0625 "<<endl
    << "--iter <n> - number of iterations, default 10 "<<endl
    << "--conductance <f> - conductance parameter, default 1.0"<<std::endl
    << "--float save image in float format"<<std::endl
    << "--short save image in short format"<<std::endl
    << "--byte save image in byte format"<<std::endl;
}


int main (int argc, char **argv)
{
  
  int clobber=0;
  int curvature=0;
  int store_float=0;
  int store_short=0;
  int store_byte=0;
  int c;
  double time_step=0.0625;
  double conductance=1.0;
  int iterations=10;
  
#if ( ITK_VERSION_MAJOR < 4 ) 
  char *_history = time_stamp(argc, argv);
  std::string history=_history;
  free(_history);
#else
  std::string history= minc_timestamp(argc,argv);
#endif

  static struct option long_options[] = { 
    {"clobber", no_argument, &clobber, 1},
    {"curvature", no_argument, &curvature, 1},
    {"gradient", no_argument, &curvature, 0},
    {"float", no_argument, &store_float, 1},
    {"short", no_argument, &store_short, 1},
    {"byte", no_argument, &store_byte, 1},
    {"step", required_argument, 0, 's'},
    {"iter", required_argument, 0, 'i'},
    {"conductance", required_argument, 0, 'c'},
    {0, 0, 0, 0}};

  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "s:i:c:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
      break;
    case 's':
      time_step = atof(optarg);
      break;
    case 'c':
      conductance = atof(optarg);
      break;
    case 'i':
      iterations = std::stoi(optarg);
      break;
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
  
  std::string input_f(argv[optind]),out_f(argv[optind+1]);
  
  if (!clobber && !access (out_f.c_str (), F_OK))
  {
    cerr << out_f.c_str () << " Exists!" << endl;
    return 1;
  }
  try
  {
#ifdef HAVE_MINC4ITK
    itk::RegisterMincIO();
#endif    
    itk::ImageFileReader<image3d >::Pointer reader = itk::ImageFileReader< image3d >::New();
    
    //initializing the reader
    reader->SetFileName(input_f.c_str());
    reader->Update();
    
    image3d::Pointer img_in=reader->GetOutput();
    
    AnisotropicFilter::Pointer filter;
    
    if(curvature)
    {
      filter = CurvatureAnisotropicFilter::New();
    } else  {
      filter = GradientAnisotropicFilter::New();
    }
    
    filter->SetNumberOfIterations(iterations);
    filter->SetTimeStep(time_step);
    filter->SetConductanceParameter(conductance);
    filter->SetInput(img_in);
    filter->Update();
    
    image3d::Pointer img_out=filter->GetOutput();
    copy_metadata(img_out,img_in);
    append_history(img_out,history);
    
#ifdef HAVE_MINC4ITK
    if(store_float)
    {
      minc::set_minc_storage_type(img_out,NC_FLOAT,true);
    } else if(store_short) {
      minc::set_minc_storage_type(img_out,NC_SHORT,true);
    } else if(store_byte) {
      minc::set_minc_storage_type(img_out,NC_BYTE,false);
    }
#else
  if(store_float)
    {
      minc::set_minc_storage_type(img_out,typeid(float).name());
    } else if(store_short) {
      minc::set_minc_storage_type(img_out,typeid(unsigned short).name());
    } else if(store_byte) {
      minc::set_minc_storage_type(img_out,typeid(unsigned char).name());
    }
#endif

    itk::ImageFileWriter< image3d >::Pointer writer = itk::ImageFileWriter<image3d >::New();
    writer->SetFileName(out_f.c_str());
    writer->SetInput( img_out );
    writer->Update();
  }
#ifdef HAVE_MINC4ITK
  catch (const generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    cerr << err.msg()<<endl;
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
