#include <iostream>

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>

#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"

typedef itk::OrientedImage<float, 3>  ImageType;
 
int main(int argc,char **argv)
{
  if(argc<2)
  {
    std::cerr<<"Usage:"<<argv[0]<<" <in.mnc> <out.mnc>"<<std::endl;
    return 1;
  }
  
  try
  {
    //registering the MINC_IO factory
    itk::ObjectFactoryBase::RegisterFactory(itk::MincImageIOFactory::New());
    /* READING */
    std::cout<<"Reading "<<argv[1]<<"..."<<std::endl;
    
    itk::ImageFileReader<ImageType >::Pointer reader = itk::ImageFileReader<ImageType >::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    
    ImageType::Pointer img=reader->GetOutput();
     
    /* WRITING */
    std::cout<<"Writing "<<argv[2]<<"..."<<std::endl;
    itk::ImageFileWriter< ImageType >::Pointer writer = itk::ImageFileWriter<ImageType >::New();
    writer->SetFileName(argv[2]);
    writer->SetInput( img );
    writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
  return 0;
}
