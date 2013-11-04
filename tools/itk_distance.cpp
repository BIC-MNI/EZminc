/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 


#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkBinaryThresholdImageFilter.h>

#ifdef HAVE_MINC4ITK
#include <time_stamp.h>    // for creating minc style history entry
#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "itkMincHelpers.h"
#else
#include "itk4MincHelpers.h"
#endif 
using namespace minc;

using namespace std;
void show_usage(const char *name)
{
  std::cerr 
    << "Usage: "<<name<<" <input> <output> " << endl
    << "--verbose be verbose "    << endl
    << "--clobber clobber output files"<<endl
    << "--signed produce signed distance map"<<endl
    << "--label <n> -extract label n before calculating distance"<<endl
    << "--voronoi <out.mnc> output voronoi map"<<std::endl;
}


typedef itk::Image<short,3> InputImageType;
typedef itk::Image<float,3> DistanceImageType;

typedef itk::BinaryThresholdImageFilter< 
                    InputImageType, 
                    InputImageType >     BinaryThresholdFilter;

typedef itk::DanielssonDistanceMapImageFilter< InputImageType, DistanceImageType >       DistanceMapFilter;
typedef itk::SignedDanielssonDistanceMapImageFilter< InputImageType, DistanceImageType > SignedDistanceMapFilter;


int main (int argc, char **argv)
{
  
  int verbose=1;
  double sigma=0.5;
  double keep=1.0;
  int order=5;
  int approx=0;
  int ss=0;
  int clobber=0;
  int label=0;
  bool label_set=false;
  std::string voronoi_f;
  
  int store_float=0;
  int store_short=0;
  int store_byte=0;
  
  
#if ( ITK_VERSION_MAJOR < 4 ) 
  char *_history = time_stamp(argc, argv);
  std::string history=_history;
  free(_history);
#else
  std::string history= minc_timestamp(argc,argv);
#endif
  
  //int voxel_neibourhood=5;
  static struct option long_options[] = { 
    {"clobber", no_argument, &clobber, 1},
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"signed",no_argument, &ss, 1},
    {"label" ,required_argument, 0, 'l'},
    {"voronoi" ,required_argument, 0, 'v'},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {0, 0, 0, 0}
  };
  
  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vq", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c)
    {
    //  case 'n':
    //    voxel_neibourhood=atoi(optarg);break;
      case 0:
        break;
      case 'l':
        label=atoi(optarg);label_set=true;break;
      case 'v':
        voronoi_f=optarg;break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
    }
  }
  if ((argc - optind) < 2) {
    show_usage (argv[0]);
    return 1;
  }
  std::string input_f=argv[optind],  out_f=argv[optind+1];
  
  // check if the file already present
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
    itk::ImageFileReader<InputImageType >::Pointer reader = itk::ImageFileReader< InputImageType >::New();
    
    //initializing the reader
    reader->SetFileName(input_f.c_str());
    reader->Update();
    
    InputImageType::Pointer input=reader->GetOutput();
    DistanceImageType::Pointer output;
#if (ITK_VERSION_MAJOR==3)
    DistanceImageType::Pointer voronoi;
#else
    SignedDistanceMapFilter::VoronoiImageType::Pointer voronoi;    
#endif
    
    BinaryThresholdFilter::Pointer threshold;
    
    if(label_set)
    {
      threshold=BinaryThresholdFilter::New();
      threshold->SetInput(input);
      threshold->SetLowerThreshold(label);
      threshold->SetUpperThreshold(label);
      threshold->Update();
      input=threshold->GetOutput();
      
    }
    
    if(ss)
    {
      SignedDistanceMapFilter::Pointer dist(SignedDistanceMapFilter::New());
      dist->SetInput(input);
      dist->Update();
      output=dist->GetOutput();
      voronoi=dist->GetVoronoiMap();
      
    } else {
      DistanceMapFilter::Pointer dist(DistanceMapFilter::New());
      dist->SetInput(input);
      dist->Update();
      output=dist->GetOutput();
      voronoi=dist->GetVoronoiMap();
    }
    

    copy_metadata(output,input);
    copy_metadata(voronoi,input);
    
    append_history(output,history);
    append_history(voronoi,history);
    
#ifdef HAVE_MINC4ITK
    if(store_float)
    {
      minc::set_minc_storage_type(output,NC_FLOAT,true);
      minc::set_minc_storage_type(voronoi,NC_FLOAT,true);
    } else if(store_short) {
      minc::set_minc_storage_type(output,NC_SHORT,true);
      minc::set_minc_storage_type(voronoi,NC_SHORT,true);
    } else if(store_byte) {
      minc::set_minc_storage_type(output,NC_BYTE,false);
      minc::set_minc_storage_type(voronoi,NC_BYTE,false);
    }
#else
  //TODO: convert to ITK4 style
#endif 
    
    itk::ImageFileWriter< DistanceImageType >::Pointer writer = itk::ImageFileWriter<DistanceImageType >::New();
    writer->SetFileName(out_f.c_str());
    writer->SetInput( output );
    writer->Update();
    
    if(!voronoi_f.empty())
    {
      
#if (ITK_VERSION_MAJOR==3)
      itk::ImageFileWriter< DistanceImageType >::Pointer writer = itk::ImageFileWriter<DistanceImageType >::New();
#else
      itk::ImageFileWriter< SignedDistanceMapFilter::VoronoiImageType >::Pointer writer = itk::ImageFileWriter<SignedDistanceMapFilter::VoronoiImageType >::New();
#endif
      writer->SetFileName(voronoi_f.c_str());
      writer->SetInput( voronoi );
      writer->Update();
    }
    
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
