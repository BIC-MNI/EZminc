#include <iostream>

#include <itkImage.h>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"

typedef itk::OrientedImage<float, 3>  ImageType;
 
int main(int argc,char **argv)
{
  if(argc<2)
  {
    std::cerr<<"Usage:"<<argv[0]<<" <in.mnc> "<<std::endl;
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
    
    // do some calculations
    ImageType::Pointer img=reader->GetOutput();
    
    const ImageType::RegionType reg=img->GetRequestedRegion();
    itk::ImageRegionConstIteratorWithIndex<ImageType> it(img, reg);
    
    double mean=0.0;
    double mx=0.0;
    double my=0.0;
    double mz=0.0;
    int cnt=0;
    double max=-1e10;
    double min=1e10;
    double max_x,max_y,max_z;
    double min_x,min_y,min_z;
    
    for(it.GoToBegin();!it.IsAtEnd();++it)
    {
      float v=it.Value();
      ImageType::IndexType idx=it.GetIndex();
      itk::Point<double,3>  p;
      img->TransformIndexToPhysicalPoint(idx,p);
      if(v>max)
      {
        max=v;
        max_x=p[0];
        max_y=p[1];
        max_z=p[2];
      }
      
      if(v<min)
      {
        min=v;
        min_x=p[0];
        min_y=p[1];
        min_z=p[2];
      }
      mean+=v;
      mx+=v*p[0];
      my+=v*p[1];
      mz+=v*p[2];
      cnt++;
    }
    if(cnt>0)
    {
      mx/=mean;
      my/=mean;
      mz/=mean;
      mean/=cnt;
    }
    
    std::cout<<"Image parameters:"<<std::endl
             <<"\tmax="<<max<<" at ("<<max_x<<","<<max_y<<","<<max_z<<")"<<std::endl
             <<"\tmin="<<min<<" at ("<<min_x<<","<<min_y<<","<<min_z<<")"<<std::endl
             <<"\tmean="<<mean<<std::endl
             <<"\tcenter of mass at ("<<mx<<","<<my<<","<<mz<<")"<<std::endl;
    
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
  return 0;
}
