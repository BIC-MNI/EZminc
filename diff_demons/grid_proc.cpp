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
#endif //HAVE_CONFIG_H

#include "mincUtils.h"
#include <itkCommand.h>

#include <itkDisplacementToVelocityFieldLogFilter.h>
#include <itkInvertDisplacementFieldImageFilter.h>
#include <itkInverseDisplacementFieldImageFilter.h>
#include <itkVectorMagnitudeImageFilter.h>
#include <itkExponentialDisplacementFieldImageFilter.h>
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <itkTransformToDisplacementFieldFilter.h>
//#include <itkDeformationFieldJacobianDeterminantFilter.h>
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include "minc/mincVectorHarmonicEnergyFilter.h"

#include <getopt.h>
#include <unistd.h>

typedef itk::Vector<float, 3 >    VectorType;
typedef itk::Image<VectorType, 3 > VectorImageType;
typedef itk::Image<float, 3 > ScalarImageType;

typedef itk::ExponentialDisplacementFieldImageFilter< VectorImageType, VectorImageType >                      _ExponentFilterType;
typedef itk::InverseDisplacementFieldImageFilter<VectorImageType, VectorImageType>                            _InvFilterType;
typedef itk::DisplacementToVelocityFieldLogFilter< VectorImageType, VectorImageType >                         _LogFilterType;
typedef itk::VectorMagnitudeImageFilter<VectorImageType,ScalarImageType>                                      _MagFilterType;
typedef itk::DisplacementFieldJacobianDeterminantFilter<VectorImageType,float>                                _DetFilterType;

typedef itk::VectorHarmonicEnergyFilter<VectorImageType,ScalarImageType>                                     _HarmonicEnergyFilterType;


typedef itk::TransformToDisplacementFieldFilter<VectorImageType,float> _Transform2DefType;


class CommandProgressUpdate : public itk::Command
{
  public:
    typedef CommandProgressUpdate   Self;
    typedef itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer;
    itkNewMacro( Self );
  protected:
    CommandProgressUpdate() {};
  public:
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      const itk::ProcessObject * filter =
          dynamic_cast< const itk::ProcessObject * >( object );
      if( ! itk::ProgressEvent().CheckEvent( &event ) )
      {
        return;
      }
      std::cout.precision(3);
      std::cout << filter->GetProgress()*100 << "%\t" <<std::flush;
    }
};


void show_usage (const char * prog)
{
  std::cerr<<"Usage:"<<prog<<" in_grid.mnc out_grid.mnc "<<std::endl
      <<"[ "<<std::endl
      <<"  --log calculate vector field logarithm"<<std::endl
      <<"  --exp calculate vector field exponential"<<std::endl
      <<"  --mag calculate vector field magnitude"<<std::endl
      <<"  --det calculate vector field Jacobian determinant"<<std::endl
      <<"  --inv calculate vector field inversion"<<std::endl
      <<"  --harm calculate harmonic energy map" << std::endl
      <<"  --iter <n> number of iterations algorithm is allowed to do"<<std::endl
      <<"  --clobber clobber output file(s)" <<std::endl
      <<"  --float save image in float format" << std::endl
      <<"  --short save image in short format" << std::endl
      <<"  --byte save image in byte format"   << std::endl
      <<"]"<<std::endl; 
}

template<class Image> void save_file(const std::string &fname,typename Image::Pointer img)
{
  typename itk::ImageFileWriter< Image >::Pointer writer = itk::ImageFileWriter<Image >::New();
  writer->SetFileName(fname.c_str());
  writer->SetInput( img );
  writer->Update();
}

template<class Image> typename Image::Pointer load_file(const std::string &fname)
{
  typename itk::ImageFileReader< Image >::Pointer reader = itk::ImageFileReader< Image >::New();
  reader->SetFileName ( fname.c_str() );
  reader->Update();
  typename Image::Pointer img = reader->GetOutput();
  img->DisconnectPipeline();
  return img;
}


int main (int argc, char **argv)
{
  int verbose=0,clobber=0;
  int glog=0,gexp=0,gsqrt=0,gsq=0,ginv=0;
  int iter=10;
  int gmag=0,gdet=0;
  int gcompat=1;
  int gharm=0;
  int store_float=0;
  int store_short=0;
  int store_byte=0;
  
  std::string minc_history=minc_timestamp(argc,argv);

  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"log",     no_argument,       &glog, 1},
    {"exp",     no_argument,       &gexp, 1},
    {"sqrt",    no_argument,       &gsqrt, 1},
    {"square",  no_argument,       &gsq, 1},
    {"inv",     no_argument,       &ginv, 1},
    {"mag",     no_argument,       &gmag, 1},
    {"det",     no_argument,       &gdet, 1},
    {"harm",    no_argument,       &gharm, 1},
    {"iter",    required_argument, 0, 'i'},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte, 1},
    {0, 0, 0, 0}
  };
  
  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vi:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c)
    {
      case 0:
        break;
      case 'i': 
        iter=atoi(optarg);
        break;
      case 'v':
        std::cout << "Version: 0.1" << std::endl;
        return 0;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
    }
  }

  if((argc - optind) < 2) {
    show_usage (argv[0]);
    return 1;
  }
  std::string out_file=argv[optind+1];
  
  if (!clobber && !access(out_file.c_str(), F_OK))
  {
    std::cerr << out_file.c_str() << " Exists!" << std::endl;
    return 1;
  }
  
  try
  {
    CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
    
    VectorImageType::Pointer  grid=load_file<VectorImageType>(argv[optind]);
    
    if(gmag)
    {
      _MagFilterType::Pointer filter=_MagFilterType::New();
      filter->SetInput(grid);
      filter->Update();
      mincify(filter->GetOutput(),minc_history,store_byte?typeid(unsigned char).name():store_short?typeid(short).name():typeid(float).name(),grid);
      save_file<ScalarImageType>(out_file,filter->GetOutput());
      
    } else if(gdet) {
      _DetFilterType::Pointer filter=_DetFilterType::New();
      filter->SetInput(grid);
      filter->Update();
      mincify(filter->GetOutput(),minc_history,store_byte?typeid(unsigned char).name():store_short?typeid(short).name():typeid(float).name(),grid);
      save_file<ScalarImageType>(out_file,filter->GetOutput());
    } else if(gharm) {
      _HarmonicEnergyFilterType::Pointer filter=_HarmonicEnergyFilterType::New();
      filter->SetInput(grid);
      filter->Update();
      mincify(filter->GetOutput(),minc_history,store_byte?typeid(unsigned char).name():store_short?typeid(short).name():typeid(float).name(),grid);
      save_file<ScalarImageType>(out_file,filter->GetOutput());
    }  else {
      VectorImageType::Pointer  out;
      if(gexp)
      {
        _ExponentFilterType::Pointer filter=_ExponentFilterType::New();
        filter->AddObserver( itk::ProgressEvent(), observer );
        filter->SetInput(grid);
        filter->Update();
        out=filter->GetOutput();
      } else if(ginv) {
        _InvFilterType::Pointer filter=_InvFilterType::New();
        filter->AddObserver( itk::ProgressEvent(), observer );
        filter->SetInput(grid);
        //filter->SetMaximumNumberOfIterations(iter);
        filter->Update();
        out=filter->GetOutput();
      } else if(glog) {
        _LogFilterType::Pointer filter=_LogFilterType::New();
        filter->AddObserver( itk::ProgressEvent(), observer );
        filter->SetInput(grid);
        //filter->SetMaximumNumberOfIterations(iter);
        filter->Update();
        out=filter->GetOutput();
      } else {
        std::cerr<<"Please specify operation!"<<std::endl;
        return 10;
      }
      
      mincify(out,minc_history,store_byte?typeid(unsigned char).name():store_short?typeid(short).name():typeid(float).name(),grid);
      save_file<VectorImageType>(out_file,out);
      std::cout << std::endl;
    }
    return 0;
  }  catch ( itk::ExceptionObject& err )
  {
    std::cout << "Could not read the fixed image information." << std::endl;
    std::cout << err << std::endl;
    return 1;
  }
  return 0;
};
