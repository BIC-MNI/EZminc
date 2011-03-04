#include <iostream>

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>

#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "minc_helpers.h"

typedef minc::def3d  ImageType;
 
int main(int argc,char **argv)
{
  if(argc<2)
  {
    std::cerr<<"Usage:"<<argv[0]<<" <in.mnc> <out.mnc>"<<std::endl;
    return 1;
  }
  
  try
  {
    /* READING */
    std::cout<<"Reading "<<argv[1]<<"..."<<std::endl;
    typedef itk::MincImageIO ImageIOType;
    ImageIOType::Pointer minc2ImageIO = ImageIOType::New();
    
    itk::ImageFileReader<ImageType >::Pointer reader = itk::ImageFileReader<ImageType >::New();
    reader->SetFileName(argv[1]);
    reader->SetImageIO( minc2ImageIO );
    reader->Update();
    
    ImageType::Pointer img=reader->GetOutput();
    std::cout<<img;
     
    /* WRITING */
    std::cout<<"Writing "<<argv[2]<<"..."<<std::endl;
    itk::ImageFileWriter< ImageType >::Pointer writer = itk::ImageFileWriter<ImageType >::New();
    writer->SetFileName(argv[2]);
    writer->SetImageIO( minc2ImageIO );
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
