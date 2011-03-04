#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "minc_helpers.h"
#include <iostream>
#include <getopt.h>
#include <math.h>

using namespace  std;
using namespace  minc;

                    
void show_usage (const char *name)
{
  std::cerr 
	  << "Usage: " << name << " <input.mnc> " << endl
    << "--verbose" << endl
    << "--mask <initial_mask.mnc> "<<endl
    << "--max-order <1|2> "<<endl;
}

int main (int argc, char **argv)
{
	int c;
  int neiborhood=10;
  int min_count=3;
  int verbose=0;
  int order=1;
  std::string out_f,mask_f;
  
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose,1},
    {"mask",  required_argument, 0, 'm'},
    {"max-order", required_argument, 0, 'o'},
    {0, 0, 0, 0}};

	for (;;)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "vm:o:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
    case 'm':
      mask_f=optarg;
      break;
    case 'o':
      order=atoi(optarg);
      break;  
		case '?':
		default:
			show_usage (argv[0]);
			return 1;
		}
	}

  if((argc - optind)<1)
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string input=argv[optind];
  
    
  try
  {
    itk::ObjectFactoryBase::RegisterFactory(itk::MincImageIOFactory::New());
    
    itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
    
    reader->SetFileName(input.c_str());
    reader->Update();
    minc::image3d::Pointer img=reader->GetOutput();
    
    mask3d::Pointer  mask;

    if(!mask_f.empty())
    {
      itk::ImageFileReader<minc::mask3d >::Pointer reader = itk::ImageFileReader<minc::mask3d >::New();
      reader->SetFileName(mask_f.c_str());
      reader->Update();
      mask=reader->GetOutput();
      
      if(!check_same(img,mask))
      {
        std::cerr<<"Error: mask & image dimensions (step size, direction cosines, dimension lengths) mismatch!"<<std::endl;
        return 1;
      }
    } else {
      mask=minc::mask3d::New();
      allocate_same(mask,img);
      mask->FillBuffer(1);
    }
    
    //save_minc(output.c_str(),mask);
    //1: find center of gravity
    minc::mask3d_iterator         it(mask, mask->GetRequestedRegion());
    minc::image3d_const_iterator  iti(img, img->GetRequestedRegion());
    
    int count=0;
    double weight=0.0;
    double wx=0.0,wy=0.0,wz=0.0;
    double voxel_volume=img->GetSpacing()[0]*img->GetSpacing()[1]*img->GetSpacing()[2];
    //minc::image3d::RegionType reg=mask->GetRequestedRegion();
    for(it.GoToBegin(),iti.GoToBegin();!it.IsAtEnd();++it,++iti)
    {
      if(!it.Value()) continue;
      tag_point pnt;
      img->TransformIndexToPhysicalPoint(iti.GetIndex(),pnt);
      float v=iti.Value()*voxel_volume;
      wx+=v*pnt[0];
      wy+=v*pnt[1];
      wz+=v*pnt[2];
      weight+=v;
      count++;
    }
    
    if(!count)
    {
      std::cout<<"0.0,0.0,0.0,0.0"<<std::endl;
      std::cerr<<"No data points (check your mask!)"<<std::endl;
      return 0;
    } else {
      wx/=weight;
      wy/=weight;
      wz/=weight;
    }
    
    std::cout<<weight<<","
             <<wx<<","
             <<wy<<","
             <<wz;
    
    if(order<=1) 
    {
      std::cout<<std::endl;
      return 0;
    }
    
    std::vector<double> moments(9,0.0);
    
    for(it.GoToBegin(),iti.GoToBegin();!it.IsAtEnd();++it,++iti)
    {
      if(!it.Value()) continue;
      tag_point pnt;
      img->TransformIndexToPhysicalPoint(iti.GetIndex(),pnt);
      float v=iti.Value()*voxel_volume/weight;
      pnt[0]-=wx;
      pnt[1]-=wy;
      pnt[2]-=wz;
      
      for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
      {
        moments[j*3+i]+=v*pnt[i]*pnt[j];
      }
    }
    for(int i=0;i<9;i++)
      std::cout<<","<<moments[i];
    
    std::cout<<std::endl;
    
		return 0;
	} catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
		cerr << err.msg()<<endl;
    return 1;
  }
	return 0;
};
