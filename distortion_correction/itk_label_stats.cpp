#include "minc_wrappers.h"
#include "minc_io.h"
#include "fftw_wrappers.h"
#include <iostream>
#include <getopt.h>
#include <math.h>
#include <limits>
#include <map>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <unistd.h>

using namespace  std;
using namespace  minc;


typedef itk::Image<unsigned char,3> Byte3DImage;
typedef itk::Image<int,3>           Int3DImage;
typedef itk::Image<float,3>         Float3DImage;

void show_usage (const char *name)
{
  std::cerr 
    << "Usage: " << name << " <input.mnc> [output]" << endl
    << "--verbose" << endl
    << "--mask <initial_mask.mnc> "<<endl
    << "--invert-mask invert mask" << endl
    << "--bg include background (label 0)"<< endl
    << "--volume <input.mnc> calculate per volume means for this volume"<<endl;
}

int main (int argc, char **argv)
{
  int c;
  int neiborhood=10;
  int min_count=3;
  int verbose=0;
  int invert_mask=0;
  int include_bg=0;
  int clobber=0;
  std::string out_f,mask_f;
  std::string vol_f;
  
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose,1},
    {"clobber", no_argument, &clobber,1},
    {"mask",    required_argument, 0, 'm'},
    {"volume",    required_argument, 0, 'v'},
    {"invert-mask",no_argument, &invert_mask,1},
    {"bg",no_argument, &include_bg,1},
    {0, 0, 0, 0}};

  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "vm:", long_options, &option_index);

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
    case 'v':
      vol_f=optarg;
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
  std::string output;
  
  if((argc - optind)>1)
  {
    output=argv[optind+1];
    if (!clobber && !access (output.c_str (), F_OK))
    {
      std::cerr << output.c_str () << " Exists!" << std::endl;
      return 1;
    }
  }
  
  try
  {
    typedef itk::ImageFileReader<Int3DImage >    ImageReaderType;
    typedef itk::ImageFileReader<Byte3DImage >   MaskReaderType;
    typedef itk::ImageFileReader<Float3DImage>   FloatImageReaderType;
    

    ImageReaderType::Pointer reader = ImageReaderType::New();

    reader->SetFileName(input.c_str());
    reader->Update();
    
    Int3DImage::Pointer img=reader->GetOutput();
    Byte3DImage::Pointer  mask(Byte3DImage::New());
    Float3DImage::Pointer vol(Float3DImage::New());
    
    if(!mask_f.empty())
    {
      MaskReaderType::Pointer reader = MaskReaderType::New();

      reader->SetFileName(mask_f.c_str());
      reader->Update();
      
      mask=reader->GetOutput();
    } else {
      allocate_same<Byte3DImage,Int3DImage>(mask,img);
      mask->FillBuffer(1-invert_mask);
    }
    
    if(!vol_f.empty())
    {
      FloatImageReaderType::Pointer reader = FloatImageReaderType::New();

      reader->SetFileName(vol_f.c_str());
      reader->Update();
      
      vol=reader->GetOutput();
    } else {
      allocate_same<Float3DImage,Int3DImage>(vol,img);
      vol->FillBuffer(0.0);
    }
    
    itk::ImageRegionConstIteratorWithIndex<Byte3DImage>   it(mask, mask->GetRequestedRegion());
    itk::ImageRegionConstIteratorWithIndex<Int3DImage>   iti(img, img->GetRequestedRegion());
    itk::ImageRegionConstIteratorWithIndex<Float3DImage> itv(vol, vol->GetRequestedRegion());
    
    int count=0;
    
    std::map<int,double> label_map;
    std::map<int,double> val_map;
    std::map<int,double> x_map,y_map,z_map;
    
    double voxel_volume=::fabs(img->GetSpacing()[0]*img->GetSpacing()[1]*img->GetSpacing()[2]);
    
    //minc::image3d::RegionType reg=mask->GetRequestedRegion();
    for(it.GoToBegin(),iti.GoToBegin(),itv.GoToBegin();!it.IsAtEnd();++it,++iti,++itv)
    {
      if(it.Value()==invert_mask) continue;
      int v=iti.Value();
      if(!include_bg && v==0) continue;
      
      tag_point pnt;
      img->TransformIndexToPhysicalPoint(iti.GetIndex(),pnt);

      label_map[v]+=voxel_volume;
      val_map[v]+=voxel_volume*itv.Value();

      x_map[v]+=pnt[0]*voxel_volume;
      y_map[v]+=pnt[1]*voxel_volume;
      z_map[v]+=pnt[2]*voxel_volume;
      count++;
    }
    
    std::ostream *out;
    
    if(!output.empty())
    {
      out = new std::ofstream(output.c_str());
    } else {
      out = &std::cout ;
    }
    
    *out<<"id,volume,mx,my,mz";
    
    if(!vol_f.empty())
      *out<<",val";
    
    *out<<std::endl;

    for(std::map<int,double>::const_iterator i=label_map.begin();i!=label_map.end();++i)
    {
      int label=(*i).first;
      double vol=(*i).second;
      *out<<label<<","
               <<vol<<","
               <<x_map[label]/vol<<","
               <<y_map[label]/vol<<","
               <<z_map[label]/vol;
               
      if(!vol_f.empty())
        *out<<val_map[label]/vol<<",";
      
      *out<<std::endl;
    }
    
    if(!output.empty())
      delete out;
    return 0;
  }
#ifdef HAVE_MINC4ITK
  catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
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
};

// kate: space-indent on; indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;show-tabs on
