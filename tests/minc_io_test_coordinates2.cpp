#include <iostream>
#include <minc_helpers.h>

typedef itk::OrientedImage<float, 3>  ImageType;
 
int main(int argc,char **argv)
{
  if(argc<2)
  {
    std::cerr<<"Usage:"<<argv[0]<<" <in.mnc> "<<std::endl;
    std::cerr<<"Will print coordinates of first 20 voxels with value above 0"<<std::endl;
    return 1;
  }
  
  try
  {
    ImageType::Pointer img=minc::load_minc<ImageType>(argv[1]);
    const ImageType::RegionType reg=img->GetRequestedRegion();
    itk::ImageRegionConstIteratorWithIndex<ImageType> it(img, reg);
    int cnt=0;
    for(it.GoToBegin();!it.IsAtEnd() && cnt<20;++it)
    {
      float v=it.Value();
      if(v>0.0)
      {
        ImageType::IndexType idx=it.GetIndex();
        itk::Point<double,3>  p;
        img->TransformIndexToPhysicalPoint(idx,p);
        std::cout<<p[0]<<","<<p[1]<<","<<p[2]<<std::endl;
        ++cnt;
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
