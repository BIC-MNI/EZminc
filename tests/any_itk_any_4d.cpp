#include <iostream>

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkNthElementImageAdaptor.h>

#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "minc_helpers.h"

typedef itk::VectorImage<float, 3>  ImageType4d;
typedef itk::OrientedImage<float, 3>  ImageType3d;
typedef itk::NthElementImageAdaptor<ImageType4d, float > ImageAdaptorType;

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
    std::cout<<"Reading "<<argv[1]<<" ..."<<std::endl;
    
    itk::ImageFileReader<ImageType4d >::Pointer reader = itk::ImageFileReader<ImageType4d >::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    
    ImageType4d::Pointer img=reader->GetOutput();
    
#if 0    
    ImageAdaptorType::Pointer adaptor=ImageAdaptorType::New();
    
    adaptor->SetImage(img);//a hack(?)
    adaptor->SelectNthElement(0);
    ImageType3d::Pointer img_o=ImageType3d::New();
    
    minc::nearest_resample_like(img_o,adaptor,adaptor,0);
    /* WRITING */
    std::cout<<"Writing "<<argv[2]<<"..."<<std::endl;
    itk::ImageFileWriter< ImageType3d >::Pointer writer = 
        itk::ImageFileWriter<ImageType3d >::New();
    
    writer->SetFileName(argv[2]);
    writer->SetInput( img_o );
    writer->Update();
#else
    std::cout<<"Writing "<<argv[2]<<" ..."<<std::endl;
    itk::ImageFileWriter< ImageType4d >::Pointer writer = 
        itk::ImageFileWriter<ImageType4d >::New();

    writer->SetFileName(argv[2]);
    writer->SetInput( img );
    writer->Update();
                   
#endif //                            
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
  return 0;
}
