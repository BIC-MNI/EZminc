// based on https://github.com/stnava/ANTs/blob/master/Utilities/itkantsReadWriteTransform.h
// with modifications from Vladimir S. FONOV 


#include <iostream>

#include <itkImage.h>
#include <itkImageIOFactory.h>
#include <itkImageIOBase.h>
#include <itkFlipImageFilter.h>
#include <itkMetaDataObject.h>
#include <itkMetaDataDictionary.h>
#include <itkCastImageFilter.h>

#include "itk4MincHelpers.h"

#include <getopt.h>
#include <unistd.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <itkDisplacementFieldTransform.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

#include <itkCompositeTransform.h>

namespace itk
{
template <class T, unsigned VImageDimension>
typename itk::Transform<T, VImageDimension, VImageDimension>::Pointer
ReadTransform(const std::string & filename) 
{
  // We must explicitly check for file existance because failed reading is an acceptable
  // state for non-displacment feilds.
  if( !itksys::SystemTools::FileExists( filename.c_str() ) )
    {
    std::cerr << "Transform file does not exist: " << filename << std::endl;
    return ITK_NULLPTR;
    }

  bool hasTransformBeenRead = false;

  typedef typename itk::DisplacementFieldTransform<T, VImageDimension> DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;
  typedef itk::ImageFileReader<DisplacementFieldType>                       DisplacementFieldReaderType;
  typename DisplacementFieldReaderType::Pointer fieldReader = DisplacementFieldReaderType::New();
  typedef typename itk::CompositeTransform<T, VImageDimension> CompositeTransformType;

  // There are known tranform type extentions that should not be considered as imaging files
  // That would be used as deformatino feilds
  // If file is an hdf5 file, assume it is a tranform instead of an image.
  if( filename.find(".h5") == std::string::npos
      && filename.find(".hdf5") == std::string::npos
      && filename.find(".hdf4") == std::string::npos
      && filename.find(".mat") == std::string::npos
      && filename.find(".txt") == std::string::npos
      && filename.find(".xfm") == std::string::npos
      )
    {
    try
      {
      fieldReader->SetFileName( filename.c_str() );
      fieldReader->Update();
      hasTransformBeenRead = true;
      }
    catch( ... )
      {
      hasTransformBeenRead = false;
      }
    }

  typedef typename itk::Transform<T, VImageDimension, VImageDimension> TransformType;
  typename TransformType::Pointer transform;
  if( hasTransformBeenRead )
    {
    typename DisplacementFieldTransformType::Pointer displacementFieldTransform =
      DisplacementFieldTransformType::New();
    displacementFieldTransform->SetDisplacementField( fieldReader->GetOutput() );
    transform = dynamic_cast<TransformType *>( displacementFieldTransform.GetPointer() );
    }
  else
    {
    typename itk::TransformFileReaderTemplate<T>::Pointer transformReader
      = itk::TransformFileReaderTemplate<T>::New();

    transformReader->SetFileName( filename.c_str() );
    try
      {
      transformReader->Update();
      }
    catch( const itk::ExceptionObject & e )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an ITK exception:\n";
      e.Print( std::cerr );
      return transform;
      }
    catch( const std::exception & e )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an exception:\n";
      std::cerr << e.what() << std::endl;
      return transform;
      }
    catch( ... )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an unknown exception!!!\n";
      return transform;
      }

    const typename itk::TransformFileReaderTemplate<T>::TransformListType * const listOfTransforms =
      transformReader->GetTransformList();
      
    if(listOfTransforms->size()>1)
    {
      typename CompositeTransformType::Pointer comp_transform=CompositeTransformType::New();
      for(typename itk::TransformFileReaderTemplate<T>::TransformListType::const_iterator i=listOfTransforms->begin();
          i!=listOfTransforms->end();
          ++i)
          {
            comp_transform->AddTransform( dynamic_cast<TransformType *>(i->GetPointer()) );
          }
      transform = dynamic_cast<TransformType *>(comp_transform.GetPointer());
    } else {
    
      transform = dynamic_cast<TransformType *>( listOfTransforms->front().GetPointer() );
    }
    
    }
  return transform;
}

template <class T, unsigned int VImageDimension>
int
WriteTransform(typename itk::Transform<T, VImageDimension, VImageDimension>::Pointer & xfrm,
               const std::string & filename)
{
  typedef typename itk::DisplacementFieldTransform<T, VImageDimension>      DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;
  typedef typename itk::ImageFileWriter<DisplacementFieldType>              DisplacementFieldWriter;
  typedef itk::TransformFileWriterTemplate<T>                               TransformWriterType;

  DisplacementFieldTransformType *dispXfrm =
    dynamic_cast<DisplacementFieldTransformType *>(xfrm.GetPointer() );

  // if it's a displacement field transform
  try
    {
    if( dispXfrm != ITK_NULLPTR )
      {
      typename DisplacementFieldType::Pointer dispField = dispXfrm->GetModifiableDisplacementField();
      typename DisplacementFieldWriter::Pointer writer = DisplacementFieldWriter::New();
      writer->SetInput(dispField);
      writer->SetFileName(filename.c_str() );
      writer->Update();
      }
    else
    // regular transform
      {
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput(xfrm);
      transformWriter->SetFileName(filename.c_str() );
      transformWriter->Update();
      }
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Can't write transform file " << filename << std::endl;
    std::cerr << "Exception Object caught: " << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace itk



void show_usage (const char *name)
{
  std::cerr  
      << "Usage: "<<name<<" <input> <output> " << std::endl
      << "\t--decompose write transform into separate files, using extension of the output file" << std::endl
      << "\t--dimensions <n> specify number of dimensions, default 3" <<std::endl
      << "\t--float - store transformation in float"<<std::endl;
      
}

int main(int argc, char ** argv)
{
  int verbose=0;
  int clobber=0;
  int dimensions=3;
  int decompose=0;
  int use_float=0;
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose, 1},
    {"quiet",   no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"decompose", no_argument, &decompose, 1},
    {"dimensions",required_argument, 0, 'd'},
    {"float",no_argument, &use_float, 1},
    {0, 0, 0, 0}
  };
    
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "d:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        break;
        
      case 'd':
        dimensions=atoi(optarg);
        break;
        
      case '?':
        /* getopt_long already printed an error message. */
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
  std::string input=argv[optind];
  std::string output=argv[optind+1];
  
  if (!clobber && !access (output.c_str(), F_OK))
  {
    std::cerr << output.c_str () << " Exists!" << std::endl;
    return 1;
  }
  
  try
  {
  
    if(use_float)
    {
      if(dimensions==3)
      {
        typedef itk::Transform<float, 3, 3> Transform;
        Transform::Pointer xfm=itk::ReadTransform<float,3>(input.c_str());
        if(xfm.IsNotNull())
        {
          if(verbose)
            std::cout<<"Tranform type="<<xfm->GetNameOfClass()<<std::endl;
          
          if(itk::WriteTransform<float,3>(xfm,output.c_str())!=EXIT_SUCCESS)
          {
            std::cerr<<"Can't write:"<<output.c_str()<<std::endl;
            return 1;
          }
        }
        else
        {
          std::cerr<<"Can't read:"<<input.c_str()<<std::endl;
          return 1;
        }
      } else if(dimensions==2) {
        typedef itk::Transform<float, 2, 2> Transform;
        
        Transform::Pointer xfm=itk::ReadTransform<float,2>(input.c_str());
        if(xfm.IsNotNull() )
        {
          if(itk::WriteTransform<float,2>(xfm,output.c_str())!=EXIT_SUCCESS)
          {
            std::cerr<<"Can't write:"<<output.c_str()<<std::endl;
            return 1;
          }
        }
        else
        {
          std::cerr<<"Can't read:"<<input.c_str()<<std::endl;
          return 1;
        }
      } else {
        std::cerr<<"Unsupported number of dimensions:"<<dimensions<<std::endl;
        return 1;
      }
    } else {
      if(dimensions==3)
      {
        typedef itk::Transform<double, 3, 3> Transform;
        Transform::Pointer xfm=itk::ReadTransform<double,3>(input.c_str());
        if(xfm.IsNotNull())
        {
          if(verbose)
            std::cout<<"Tranform type="<<xfm->GetNameOfClass()<<std::endl;
          
          if(itk::WriteTransform<double,3>(xfm,output.c_str())!=EXIT_SUCCESS)
          {
            std::cerr<<"Can't write:"<<output.c_str()<<std::endl;
            return 1;
          }
        }
        else
        {
          std::cerr<<"Can't read:"<<input.c_str()<<std::endl;
          return 1;
        }
      } else if(dimensions==2) {
        typedef itk::Transform<double, 2, 2> Transform;
        
        Transform::Pointer xfm=itk::ReadTransform<double,2>(input.c_str());
        if(xfm.IsNotNull() )
        {
          if(itk::WriteTransform<double,2>(xfm,output.c_str())!=EXIT_SUCCESS)
          {
            std::cerr<<"Can't write:"<<output.c_str()<<std::endl;
            return 1;
          }
        }
        else
        {
          std::cerr<<"Can't read:"<<input.c_str()<<std::endl;
          return 1;
        }
      } else {
        std::cerr<<"Unsupported number of dimensions:"<<dimensions<<std::endl;
        return 1;
      }
    }    
  }

  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }  
  
  return 0;
}