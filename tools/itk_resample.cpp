/* ----------------------------- MNI Header -----------------------------------
@NAME       :  itk_resample
@DESCRIPTION:  an example of using spline itnerpolation with MINC xfm transform
@COPYRIGHT  :
              Copyright 2011 Vladimir Fonov, McConnell Brain Imaging Centre, 
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

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>
#include <algorithm>
#include <stdlib.h>

#include <itkVector.h>
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include <itkIdentityTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkVectorImage.h>
#include <itkBinaryThresholdImageFilter.h>

#if (ITK_VERSION_MAJOR<5)
#include <itkVectorResampleImageFilter.h>
#endif

#include "mincVectorBSplineInterpolate.h"
#include <itkImageConstIterator.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkAntiAliasBinaryImageFilter.h>

#include <unistd.h>
#include <getopt.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkImageIOBase.h>

#if ( ITK_VERSION_MAJOR < 4 )
#include <vxl/core/vnl/vnl_cross.h>
#else
#include <vnl/vnl_cross.h>
#include <itkCompositeTransform.h>
#endif

#ifdef HAVE_MINC4ITK
#include <time_stamp.h>    // for creating minc style history entry
#include "itkMincHelpers.h"
#include "itkMincGeneralTransform.h"
#include "itkMincImageIO.h"
#include "itkMincHelpers.h"
#else
#include "itk4MincHelpers.h"
#include "itkMINCTransformAdapter.h"
#endif //HAVE_MINC4ITK

typedef itk::ImageBase<3>           Image3DBase;
typedef itk::Image<float,3>         Float3DImage;
typedef itk::Image<int,3>           Int3DImage;
typedef itk::Image<short,3>         Short3DImage;
typedef itk::Image<unsigned char,3> Byte3DImage;
typedef itk::Image<itk::Vector<float,3>, 3>  Vector3DImage;


typedef itk::ImageIOBase          IOBase;
typedef itk::SmartPointer<IOBase> IOBasePointer;

typedef itk::BSplineInterpolateImageFunction< Float3DImage, double, double >  InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction< Int3DImage, double >    NNInterpolatorType;
//typedef minc::mincVectorBSplineInterpolate<Vector3DImage,double>            VectorInterpolatorType;
typedef itk::LinearInterpolateImageFunction<Vector3DImage, double>            VectorInterpolatorType; // TODO: fix this for ITKv5
//typedef itk::BSplineInterpolateImageFunction< Vector3DImage, double, double >  VectorInterpolatorType;


typedef itk::ResampleImageFilter<Float3DImage, Float3DImage> FloatFilterType;
typedef itk::ResampleImageFilter<Int3DImage  , Int3DImage>   IntFilterType;

// #if ( ITK_VERSION_MAJOR < 5 )
// typedef itk::VectorResampleImageFilter<Vector3DImage , Vector3DImage>  VectorFilterType;
// #else
// typedef itk::ResampleImageFilter<Vector3DImage , Vector3DImage>  VectorFilterType;
// #endif

typedef itk::Transform<double,3,3>           TransformBaseType;
typedef itk::IdentityTransform<double,3>     IdentityTransformType;
  
#if ( ITK_VERSION_MAJOR < 4 )
typedef minc::XfmTransform<double,3,3>  TransformType;
#else
typedef itk::MINCTransformAdapter<double,3,3> TransformType;
#endif

using namespace  std;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" <input> <output.mnc> " << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--like <example> (default behaviour analogous to use_input_sampling)"<<std::endl
    << "--transform <xfm_transform> "<<std::endl
    << "--order <n> spline order, default 2 "<<std::endl
    << "--uniformize <step> - will make a volume with uniform step size and no direction cosines" << std::endl
    << "--unistep <step> - will make a volume with the same step size and preserve direction cosines" << std::endl
    << "--normalize - outputs volume with right-hand direction cosines and positive step-sizes, uses NN interpolation" << std::endl
    << "--transform-bbox apply transformation to the bounding box for the final sampling, use with uniformize" << std::endl
    << "--invert_transform  - apply inverted transform"<<std::endl
    << "--labels - assume that input data is discrete labels, will use Nearest-Neighbour interpolator"<<std::endl
    << "--byte  - store image in byte  voxels minc file"<<std::endl
    << "--short - store image in short voxels minc file"<<std::endl
    << "--float - store image in float voxels minc file"<<std::endl
    << "--relabel map.txt apply relabeling map"<<std::endl
    << "--lut map.txt  apply relabeling map"<<std::endl
    << "--lut-string \"a b;c d;....\" apply lut string in command line"<<std::endl
    << "--aa <fwhm> apply anti-aliasing filter to labels before resampling  (only usable for order > 0 )"<<std::endl
    << "--baa apply binary anti-aliasing filter to labels before resampling (only usable for order > 0 )"<<std::endl
    << "--fill <v> - fill value for voxels falling outside of input volume, default 0"<<std::endl;
}


struct resample_options
{
  int    invert;
  double uniformize;
  double unistep;
  int    transform_bbox;
  
  int    normalize;
  
  std::string history;
  int store_float;
  int store_short;
  int store_byte;
  std::map<int,int> label_map;
  
  double aa_fwhm;
  int baa_smooth;
  double fill;
  
  resample_options():
    invert(0),
    uniformize(0),unistep(0),
    transform_bbox(0),
    store_float(0),store_short(0),store_byte(0),
    aa_fwhm(0), baa_smooth(0) ,normalize(0),fill(0)
  {}
};

template<class T,class I> void generate_uniform_sampling(T* flt, const I* img, double step, TransformBaseType * tfm)
{
  //obtain physical coordinats of all image corners
  typename I::RegionType r=img->GetLargestPossibleRegion();
  std::vector<double> corner[3];

  if( step==0.0 ) //assume we want minimum reasonable step size
    step=std::min<double>(std::min<double>(
                          ::fabs(img->GetSpacing()[0]),
                          ::fabs(img->GetSpacing()[1]) ),
                          ::fabs(img->GetSpacing()[2]) );
    
  //walk over all corners
  for(int i=0; i<8; i++)
  {
    typename I::IndexType idx;
    typename I::PointType c,o;
    
    idx[0]=r.GetIndex()[0]+r.GetSize()[0]*(  i    % 2 );
    idx[1]=r.GetIndex()[1]+r.GetSize()[1]*( (i/2) % 2 );
    idx[2]=r.GetIndex()[2]+r.GetSize()[2]*( (i/4) % 2 );
    
    img->TransformIndexToPhysicalPoint(idx, c );
    
    o=tfm->TransformPoint( c );
    
    for(int j=0; j<3; j++)
      corner[j].push_back( o[j] );
  }
  
  typename I::IndexType start;
  typename T::SizeType size;
  typename T::OriginPointType org;
  typename I::SpacingType spc;
  spc.Fill(step);
  
  for(int j=0;j<3;j++)
  {
    std::sort(corner[j].begin(), corner[j].end());
    
    size[j]=ceil( (corner[j][7]-corner[j][0])/step );
    
    org[j]=corner[j][0];
  }
  
  Float3DImage::DirectionType identity;
  identity.SetIdentity();
  flt->SetOutputDirection(identity);
  start.Fill(0);
  
  flt->SetOutputStartIndex( start );
  flt->SetSize( size );
  flt->SetOutputOrigin( org );
  flt->SetOutputSpacing( spc );
}

template<class T,class I> void generate_unistep_sampling(T* flt, const I* img,double step)
{
  //obtain physical coordinats of all image corners
  typename I::RegionType r=img->GetLargestPossibleRegion();
  
  typename I::IndexType start;
  typename T::SizeType size;
  typename T::OriginPointType org=img->GetOrigin();
  typename I::SpacingType spc;
  spc.Fill(step);
  
  for(int j=0;j<3;j++)
  {
    org[j]=org[j]-img->GetSpacing()[j]/2.0+step/2.0;
    size[j]=::ceil((r.GetSize()[j]+1)*img->GetSpacing()[j]/step);
  }
  
  start.Fill(0);
  
  flt->SetOutputDirection(img->GetDirection());
  
  flt->SetOutputStartIndex(start);
  flt->SetSize(size);
  flt->SetOutputOrigin(org);
  flt->SetOutputSpacing(spc);
}


template<class D> bool abs_compare(const D& a,const D&b) { return fabs(a)<fabs(b);}

template<class T,class I> void generate_normalized_sampling(T* flt, const I* img)
{
  //obtain physical coordinats of all image corners
  typedef typename I::PointType  PointType;
  typedef typename I::IndexType  IndexType;
  typedef typename I::PointType::VectorType VectorType;
  typename I::RegionType r=img->GetLargestPossibleRegion();
  
  int i[3];
  //check all combinations of corners
  
  for(i[0]=0;i[0]<2;i[0]++)
    for(i[1]=0;i[1]<2;i[1]++)
      for(i[2]=0;i[2]<2;i[2]++)
  {
    IndexType idx,idx2[3];
    PointType corner,c2[3];
    VectorType v[3];
    double   dist[3];
    int      remap[3],rrmap[3];
    double   proj[3];
    
    for(int j=0;j<3;j++)
    {
      idx[j]=r.GetIndex()[j]+r.GetSize()[j]*i[j];
      remap[j]=-1;
      rrmap[j]=-1;
      proj[j]=0.0;
    }
    img->TransformIndexToPhysicalPoint(idx,corner);
    
    for(int j=0;j<3;j++)
    {
      idx2[j]=idx;
      idx2[j][j]=r.GetIndex()[j]+r.GetSize()[j]*(i[j]^1);
      
      img->TransformIndexToPhysicalPoint(idx2[j],c2[j]);
      
      dist[j]=corner.EuclideanDistanceTo(c2[j]);
      
      v[j]=(c2[j]-corner)/dist[j];
      
      //sort according to orientation
      size_t mp=std::max_element(v[j].GetDataPointer(),v[j].GetDataPointer()+3,
                                 abs_compare<typename VectorType::RealValueType>)-(v[j].GetDataPointer());
      
      if(rrmap[mp]==-1) //we have already found this axis?
      {
        remap[j]=mp;
        rrmap[mp]=j;
        proj[j]=v[j][mp];
      } else {
        std::cerr<<"Failed to find max component!"<<std::endl;
        std::cerr<<"mp="<<mp<<std::endl;
      }
    }
    //check if all axis are mapped
    for(int j=0;j<3;j++)
    {
      if(remap[j]==-1)//assign next available axis
      {
        int k;
        for(k=0; k<3&&rrmap[k]!=-1 ;k++);
        if(k==3) //abort?
          throw itk::ExceptionObject("Wierd orientation, can't figure out how to normalize it");
        remap[j]=k;
        rrmap[k]=j;
      }
    }
    
    vnl_vector_fixed<typename VectorType::RealValueType,3> xvec=v[rrmap[0]].GetVnlVector();
    vnl_vector_fixed<typename VectorType::RealValueType,3> yvec=v[rrmap[1]].GetVnlVector();
    vnl_vector_fixed<typename VectorType::RealValueType,3> zvec=v[rrmap[2]].GetVnlVector();
    
    //check if coordinate system is positive (right handed) i.e determinant is positive
    //and that projections on axes are posive
    if(dot_product( vnl_cross_3d(xvec,yvec), zvec ) > 0 &&
       proj[0]>0 &&
       proj[1]>0 &&
       proj[2]>0 )
    {
      //found best candidate
      //now calculate the bounds
      
      typename I::IndexType start;
      typename T::SizeType size;
      typename I::SpacingType spc;
      
      Float3DImage::DirectionType dir;
      for(int j=0;j<3;j++)
      {
        dir[j][0]=xvec[j];
        dir[j][1]=yvec[j];
        dir[j][2]=zvec[j];
        
        size[j]=r.GetSize()[rrmap[j]];
        
        spc[j]=fabs( img->GetSpacing()[rrmap[j]] ); // have to make steps positive!
        
      }
      flt->SetOutputDirection(dir);
      
      start.Fill(0);
      flt->SetOutputStartIndex(start);
      
      flt->SetSize(size);
      flt->SetOutputOrigin(corner);
      flt->SetOutputSpacing(spc);
      
      return;
    }
  }
  // should report that didn't find anything ?
}



template<class Image,class ImageOut,class Interpolator> 
void resample_image(
   IOBase* base,  
   const std::string& output_f,
   const std::string& xfm_f,
   const std::string& like_f,
   Interpolator* interpolator,
   const resample_options &opt
  )
{
  typedef typename itk::ResampleImageFilter<Image, ImageOut> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >              ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >           ImageWriterType;
  typedef typename Image::PixelType                          InputPixelType;
  typename ImageReaderType::Pointer                 reader = ImageReaderType::New();
  
  //initializing the reader
  reader->SetImageIO(base);
  reader->SetFileName(base->GetFileName());
  reader->Update();
  
  typename Image::Pointer in=reader->GetOutput();
  
  if(!opt.label_map.empty())
  {
    for(itk::ImageRegionIterator<Image> it(in, in->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
      std::map<int,int>::const_iterator pos=opt.label_map.find(static_cast<int>(it.Get()));
      if(pos==opt.label_map.end())
        it.Set(0); //set to BG
      else
        it.Set(static_cast<InputPixelType>((*pos).second));
    }
  }

  typename ResampleFilterType::Pointer filter  = ResampleFilterType::New();
  
  //creating coordinate transformation objects
  TransformType::Pointer minc_transform = TransformType::New();
  IdentityTransformType::Pointer identity_transform = IdentityTransformType::New();
  TransformBaseType* transform;
  
  if(!xfm_f.empty())
  {
#ifdef HAVE_MINC4ITK
    //reading a minc style xfm file
    minc_transform->OpenXfm(xfm_f.c_str());
#else
    minc_transform->OpenXfm(xfm_f.c_str());
#endif
    transform=minc_transform.GetPointer();
  } else {
    filter->SetTransform( identity_transform );
    transform=identity_transform.GetPointer();
  }

  //creating the interpolator
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( opt.fill );
  //this is for processing using batch system
  filter->SetNumberOfThreads(1);

  typename Image::Pointer like=nullptr;

  if( opt.transform_bbox && !xfm_f.empty() && opt.invert)
    minc_transform->Invert();
  
  if(!like_f.empty())
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();
    
    if(opt.uniformize!=0.0 || opt.transform_bbox )
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,
        reader->GetOutput(),opt.uniformize,
        opt.transform_bbox ? transform : (TransformBaseType*)identity_transform.GetPointer() );
      
    } else if(opt.unistep!=0.0) {
      generate_unistep_sampling<ResampleFilterType,Image>(filter,reader->GetOutput(),opt.unistep);
    } else if(opt.normalize) {
      generate_normalized_sampling<ResampleFilterType,Image>(filter,reader->GetOutput());
    } else {
      filter->SetOutputParametersFromImage(reader->GetOutput());
      filter->SetOutputDirection(reader->GetOutput()->GetDirection());
    }
    like=reader->GetOutput();
    like->DisconnectPipeline();
  } else {
    if(opt.uniformize!=0.0 || opt.transform_bbox)
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,
        in, opt.uniformize,
        opt.transform_bbox ? transform : (TransformBaseType*)identity_transform.GetPointer() );
    } else if(opt.unistep!=0.0) {
      generate_unistep_sampling<ResampleFilterType,Image>(filter,in,opt.unistep);
    } else if(opt.normalize) {
      generate_normalized_sampling<ResampleFilterType,Image>(filter,in);
    } else {
      //we are using original sampling
      filter->SetOutputParametersFromImage(in);
      filter->SetOutputDirection(in->GetDirection());
    }
  }
  
  if( opt.transform_bbox && !xfm_f.empty() && opt.invert)
    minc_transform->Invert();
  
  if(!opt.invert && !xfm_f.empty() ) minc_transform->Invert(); //should be inverted by default to walk through target space

  filter->SetTransform( transform );
  
  filter->SetInput(in);

#ifdef _DEBUG  
  std::cout<<"Filter:"<<filter<<std::endl;
  
  std::cout<<"Transform:"<<transform<<std::endl;

  std::cout<<"Interpolator:"<< interpolator<<std::endl;
  
  std::cout<<"Input:"<< in <<std::endl;
#endif
  
  filter->Update();
  typename ImageOut::Pointer out=filter->GetOutput();
  
  minc::copy_metadata(out,in);
  minc::append_history(out,opt.history);
  
  //correct dimension order
  if(like.IsNotNull() && opt.uniformize==0.0 && !opt.transform_bbox && !opt.normalize && opt.unistep==0.0)
    minc::copy_dimorder(out,like);
  else if( opt.uniformize!=0.0 || opt.transform_bbox || opt.normalize || opt.unistep!=0.0)
    minc::delete_dimorder(out);
    
  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  
  if( getenv("MINC_COMPRESS") != NULL)
    writer->SetUseCompression( true );
  
#ifdef HAVE_MINC4ITK
  if(opt.store_float)
  {
    minc::set_minc_storage_type(out,NC_FLOAT,true);
  } else if(opt.store_short) {or
    minc::set_minc_storage_type(out,NC_SHORT,true);
  } else if(opt.store_byte) {
    minc::set_minc_storage_type(out,NC_BYTE,false);
  }
#else
  if(opt.store_float)
  {
    minc::set_minc_storage_type(out,typeid(float).name());
  } else if(opt.store_short) {
    minc::set_minc_storage_type(out,typeid(unsigned short).name());
  } else if(opt.store_byte) {
    minc::set_minc_storage_type(out,typeid(unsigned char).name());
  }
#endif

  writer->SetInput( out );
  writer->Update();
}

template<class Image,class TmpImage,class ImageOut,class Interpolator> 
void resample_label_image (
   IOBase* base,
   const std::string& output_f,
   const std::string& xfm_f,
   const std::string& like_f,
   Interpolator* interpolator,
   const resample_options &opt )
{
  typedef typename itk::ResampleImageFilter<TmpImage, TmpImage> ResampleFilterType;
  typedef typename itk::ImageFileReader<Image >                 ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >              ImageWriterType;
  typedef itk::ImageRegionConstIterator<Image>                  ConstInputImageIteratorType;
  typedef itk::ImageRegionConstIterator<TmpImage>               ConstTmpImageIteratorType;
  typedef itk::ImageRegionIterator<TmpImage>                    TmpImageIteratorType;
  
  typedef itk::ImageRegionIterator<ImageOut>                    ImageOutIteratorType;
  typedef itk::BinaryThresholdImageFilter<Image,TmpImage>       ThresholdFilterType;
  typedef itk::SmoothingRecursiveGaussianImageFilter<TmpImage,TmpImage>  BlurFilterType;
  typedef itk::AntiAliasBinaryImageFilter<TmpImage,TmpImage>    BAAFilterType;

  typedef typename Image::PixelType InputPixelType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();

  double aa_sigma=opt.aa_fwhm/(2.0*sqrt(2*log(2.0)));
  
  //initializing the reader
  reader->SetImageIO(base);
  reader->SetFileName(base->GetFileName());
  reader->Update();

  typename Image::Pointer in=reader->GetOutput();
  
  if(!opt.label_map.empty())
  {
    for(itk::ImageRegionIterator<Image> it(in, in->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
      std::map<int,int>::const_iterator pos=opt.label_map.find(static_cast<int>(it.Get()));
      if(pos==opt.label_map.end())
        it.Set(0); //set to BG
      else
        it.Set(static_cast<InputPixelType>((*pos).second));
    }
  }

  typename ResampleFilterType::Pointer  filter           = ResampleFilterType::New();
  typename ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
  typename BlurFilterType::Pointer      blur_filter     = BlurFilterType::New();
  typename BAAFilterType::Pointer       baa_filter      = BAAFilterType::New();

  //creating coordinate transformation objects
  TransformType::Pointer minc_transform = TransformType::New();
  IdentityTransformType::Pointer identity_transform = IdentityTransformType::New();
  TransformBaseType* transform;

  if(!xfm_f.empty())
  {
#ifdef HAVE_MINC4ITK
    //reading a minc style xfm file
    minc_transform->OpenXfm(xfm_f.c_str());
#else
    minc_transform->OpenXfm(xfm_f.c_str());
#endif
    transform=minc_transform.GetPointer();
  } else {
    filter->SetTransform( identity_transform );
    transform=identity_transform.GetPointer();
  }

  //creating the interpolator
  filter->SetInterpolator( interpolator );
  filter->SetDefaultPixelValue( -9 );
  //this is for processing using batch system
  filter->SetNumberOfThreads(1);

  typename Image::Pointer like=nullptr;

  if( opt.transform_bbox && !xfm_f.empty() && opt.invert)
    minc_transform->Invert();

  if(!like_f.empty())
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();

    if(opt.uniformize!=0.0 || opt.transform_bbox )
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,
        reader->GetOutput(),opt.uniformize,
        opt.transform_bbox ? transform : (TransformBaseType*)identity_transform.GetPointer() );

    } else if(opt.unistep!=0.0) {
      generate_unistep_sampling<ResampleFilterType,Image>(filter,reader->GetOutput(),opt.unistep);
    } else if(opt.normalize) {
      generate_normalized_sampling<ResampleFilterType,Image>(filter,reader->GetOutput());
    } else {
      filter->SetOutputParametersFromImage(reader->GetOutput());
      filter->SetOutputDirection(reader->GetOutput()->GetDirection());
    }
    like=reader->GetOutput();
    like->DisconnectPipeline();
  } else {
    if(opt.uniformize!=0.0 || opt.transform_bbox)
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,
        in, opt.uniformize,
        opt.transform_bbox ? transform : (TransformBaseType*)identity_transform.GetPointer() );
    } else if(opt.unistep!=0.0) {
      generate_unistep_sampling<ResampleFilterType,Image>(filter,in,opt.unistep);
    } else if(opt.normalize) {
      generate_normalized_sampling<ResampleFilterType,Image>(filter,in);
    } else {
      //we are using original sampling
      filter->SetOutputParametersFromImage(in);
      filter->SetOutputDirection(in->GetDirection());
    }
  }

  if( opt.transform_bbox && !xfm_f.empty() && opt.invert)
    minc_transform->Invert();

  if(!opt.invert && !xfm_f.empty() ) minc_transform->Invert(); //should be inverted by default to walk through target space

  filter->SetTransform( transform );
  
  typename TmpImage::RegionType region;
  region.SetSize (filter->GetSize());
  region.SetIndex(filter->GetOutputStartIndex());
  
  
  typename TmpImage::Pointer MaxImage= TmpImage::New();
  MaxImage->SetOrigin(filter->GetOutputOrigin());
  MaxImage->SetSpacing(filter->GetOutputSpacing());
  MaxImage->SetDirection(filter->GetOutputDirection());
  
  MaxImage->SetLargestPossibleRegion(region);
  MaxImage->SetBufferedRegion(region);
  MaxImage->SetRequestedRegion(region);
  MaxImage->Allocate();
  MaxImage->FillBuffer(-10.0);
  
  typename ImageOut::Pointer LabelImage= ImageOut::New();
  LabelImage->SetOrigin(filter->GetOutputOrigin());
  LabelImage->SetSpacing(filter->GetOutputSpacing());
  LabelImage->SetDirection(filter->GetOutputDirection());
  
  LabelImage->SetLargestPossibleRegion(region);
  LabelImage->SetBufferedRegion(region);
  LabelImage->SetRequestedRegion(region);
  LabelImage->Allocate();
  
  //split the input image
  std::set<InputPixelType> sval;
  for(ConstInputImageIteratorType it(in, in->GetBufferedRegion()); !it.IsAtEnd(); ++it)
  {
    InputPixelType val = it.Get();

    if(isfinite(val))
      sval.insert(val);
  }
  
  // use default label 
  LabelImage->FillBuffer(*sval.begin());
  
  //std::vector<typename TmpImage::Pointer> list_of_labels;
  // Generate a bunch of binary copies (unfortunately, we store them as double)
  for(typename set<InputPixelType>::iterator label_it = sval.begin(); label_it != sval.end(); ++label_it)
  {
    threshold_filter->SetLowerThreshold(*label_it);
    threshold_filter->SetUpperThreshold(*label_it);
    threshold_filter->SetInsideValue(1.0);
    threshold_filter->SetOutsideValue(0.0);
    threshold_filter->SetInput(in);
    
    if(opt.aa_fwhm>0.0)
    {
      //blur_filter->SetUseImageSpacing(true);
      blur_filter->SetSigma(aa_sigma);
      blur_filter->SetInput(threshold_filter->GetOutput());
      filter->SetInput(blur_filter->GetOutput());
    } else if(opt.baa_smooth) {
      //baa_filter->SetMaximumRMSChange(baa_smooth);
      baa_filter->SetInput(threshold_filter->GetOutput());
      filter->SetInput(baa_filter->GetOutput());
    } else {
      filter->SetInput(threshold_filter->GetOutput());
    }
    
    filter->Update();
    
    ConstTmpImageIteratorType it_filter(filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
    TmpImageIteratorType      it_max(MaxImage, MaxImage->GetBufferedRegion());
    ImageOutIteratorType      it_out(LabelImage, LabelImage->GetBufferedRegion());
    
    for(; !it_filter.IsAtEnd(); ++it_filter,++it_max,++it_out)
    {
      typename TmpImage::PixelType val = it_filter.Get();
      typename TmpImage::PixelType val_max = it_max.Get();
      
      if( !isfinite(val) || val<-8.5) {
        it_max.Set(-8.5);
        it_out.Set(opt.fill);
      } else if( label_it == sval.begin() || isfinite(val) && val>val_max)
      {
        it_max.Set(val);
        it_out.Set(*label_it);
      }
    }
  }
  
  
  //typename ImageOut::Pointer out=filter->GetOutput();
  minc::copy_metadata(LabelImage,in);
  minc::append_history(LabelImage,opt.history);
  
  //correct dimension order
  if(like.IsNotNull() && opt.uniformize==0.0 && !opt.transform_bbox && !opt.normalize && opt.unistep==0.0)
    minc::copy_dimorder(LabelImage,like);
  else if( opt.uniformize!=0.0 || opt.transform_bbox || opt.normalize || opt.unistep!=0.0)
    minc::delete_dimorder(LabelImage);
    
  //generic file writer
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  
#ifdef HAVE_MINC4ITK
  if(opt.store_float)
  {
    minc::set_minc_storage_type(LabelImage,NC_FLOAT,true);
  } else if(opt.store_short) {
    minc::set_minc_storage_type(LabelImage,NC_SHORT,true);
  } else if(opt.store_byte) {
    minc::set_minc_storage_type(LabelImage,NC_BYTE,false);
  }
#else  
  if(opt.store_float)
  {
    minc::set_minc_storage_type(LabelImage,typeid(float).name());
  } else if(opt.store_short) {
    minc::set_minc_storage_type(LabelImage,typeid(unsigned short).name());
  } else if(opt.store_byte) {
    minc::set_minc_storage_type(LabelImage,typeid(unsigned char).name());
  }
#endif

  writer->SetInput( LabelImage );
  
  if( getenv("MINC_COMPRESS") != NULL)
    writer->SetUseCompression( true );
  
  writer->Update();
}



template<class Image,class ImageOut,class Interpolator> 
void resample_vector_image(
   IOBase* base,  
   const std::string& output_f,
   const std::string& xfm_f,
   const std::string& like_f,
   Interpolator* interpolator,
   const resample_options &opt
   )
{
  #if (ITK_VERSION_MAJOR<5)
  typedef typename itk::VectorResampleImageFilter<Image, ImageOut> ResampleFilterType;
  #else
  typedef typename itk::ResampleImageFilter<Image, ImageOut> ResampleFilterType;
  #endif

  typedef typename itk::ImageFileReader<Image >   ImageReaderType;
  typedef typename itk::ImageFileWriter<ImageOut >   ImageWriterType;
  
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  
  //initializing the reader
  reader->SetImageIO(base);
  reader->SetFileName(base->GetFileName());
  reader->Update();
  
  typename Image::Pointer in=reader->GetOutput();

  typename ResampleFilterType::Pointer filter  = ResampleFilterType::New();
  
  //creating coordinate transformation objects
  TransformType::Pointer minc_transform = TransformType::New();
  IdentityTransformType::Pointer identity_transform = IdentityTransformType::New();
  TransformBaseType* transform;

  if(!xfm_f.empty())
  {
#ifdef HAVE_MINC4ITK
    //reading a minc style xfm file
    minc_transform->OpenXfm(xfm_f.c_str());
#else
    minc_transform->OpenXfm(xfm_f.c_str());
#endif
    transform=minc_transform.GetPointer();
  } else {
    filter->SetTransform( identity_transform );
    transform=identity_transform.GetPointer();
  }

  //creating the interpolator
  filter->SetInterpolator( interpolator );
  //filter->SetDefaultPixelValue( 0 );
  //this is for processing using batch system
  filter->SetNumberOfThreads(1);

  typename Image::Pointer like=nullptr;

  if( opt.transform_bbox && !xfm_f.empty() && opt.invert)
    minc_transform->Invert();

  if(!like_f.empty())
  {
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(like_f.c_str());
    reader->Update();

    if(opt.uniformize!=0.0 || opt.transform_bbox )
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,
        reader->GetOutput(),opt.uniformize,
        opt.transform_bbox ? transform : (TransformBaseType*)identity_transform.GetPointer() );

    } else if(opt.unistep!=0.0) {
      generate_unistep_sampling<ResampleFilterType,Image>(filter,reader->GetOutput(),opt.unistep);
    } else if(opt.normalize) {
      generate_normalized_sampling<ResampleFilterType,Image>(filter,reader->GetOutput());
    } else {
      filter->SetOutputSpacing(reader->GetOutput()->GetSpacing());
      filter->SetOutputOrigin(reader->GetOutput()->GetOrigin());
      filter->SetOutputDirection(reader->GetOutput()->GetDirection());
    }
    like=reader->GetOutput();
    like->DisconnectPipeline();
  } else {
    if(opt.uniformize!=0.0 || opt.transform_bbox)
    {
      generate_uniform_sampling<ResampleFilterType,Image>(filter,
        in, opt.uniformize,
        opt.transform_bbox ? transform : (TransformBaseType*)identity_transform.GetPointer() );
    } else if(opt.unistep!=0.0) {
      generate_unistep_sampling<ResampleFilterType,Image>(filter,in,opt.unistep);
    } else if(opt.normalize) {
      generate_normalized_sampling<ResampleFilterType,Image>(filter,in);
    } else {
      //we are using original sampling
      filter->SetOutputSpacing(in->GetSpacing());
      filter->SetOutputOrigin(in->GetOrigin());
      filter->SetOutputDirection(in->GetDirection());
    }
  }

  if( opt.transform_bbox && !xfm_f.empty() && opt.invert)
    minc_transform->Invert();

  if(!opt.invert && !xfm_f.empty() ) minc_transform->Invert(); //should be inverted by default to walk through target space

  filter->SetTransform( transform );
  
  filter->SetInput(in);
  filter->Update();
  typename ImageOut::Pointer out=filter->GetOutput();
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(output_f.c_str());
  minc::copy_metadata(out,in);
  minc::append_history(out,opt.history);
  
  //correct dimension order
  if(like.IsNotNull() && opt.uniformize==0.0 && !opt.transform_bbox && !opt.normalize && opt.unistep==0.0)
    minc::copy_dimorder(out,like);
  else if( opt.uniformize!=0.0 || opt.transform_bbox || opt.normalize || opt.unistep!=0.0)
    minc::delete_dimorder(out);
  
  //generic file writer
  
#ifdef HAVE_MINC4ITK
  if(opt.store_float)
  {
    minc::set_minc_storage_type(out,NC_FLOAT,true);
  } else if(opt.store_short) {
    minc::set_minc_storage_type(out,NC_SHORT,true);
  } else if(opt.store_byte) {
    minc::set_minc_storage_type(out,NC_BYTE,false);
  }
#else
  if(opt.store_float)
  {
    minc::set_minc_storage_type(out,typeid(float).name());
  } else if(opt.store_short) {
    minc::set_minc_storage_type(out,typeid(unsigned short).name());
  } else if(opt.store_byte) {
    minc::set_minc_storage_type(out,typeid(unsigned char).name());
  }
#endif

  writer->SetInput( out );
  
  if( getenv("MINC_COMPRESS") != NULL)
    writer->SetUseCompression( true );
  
  writer->Update();
}


int main (int argc, char **argv)
{
  
  
  resample_options opt;
    
  int verbose=0, clobber=0;
  int order=2;
  int labels=0;
  
  std::string like_f,xfm_f,output_f,input_f;
  std::string map_f,map_str;
  
#ifdef HAVE_MINC4ITK
  char *_history = time_stamp(argc, argv); 
  opt.history=_history;
  free(_history);
#else
  opt.history = minc_timestamp(argc,argv);
#endif  
  bool order_was_set=false;
  
  static struct option long_options[] = {
    {"verbose",      no_argument,       &verbose, 1},
    {"quiet",        no_argument,       &verbose, 0},
    {"clobber",      no_argument,       &clobber, 1},
    {"like",         required_argument, 0, 'l'},
    {"transform",    required_argument, 0, 't'},
    {"order",        required_argument, 0, 'o'},
    {"uniformize",   required_argument, 0, 'u'},
    {"unistep",      required_argument, 0, 'U'},
    {"invert_transform", no_argument,   &opt.invert, 1},
    {"normalize",    no_argument,       &opt.normalize, 1},
    {"labels",       no_argument,       &labels, 1},
    {"float",        no_argument,       &opt.store_float, 1},
    {"short",        no_argument,       &opt.store_short, 1},
    {"byte",         no_argument,       &opt.store_byte,  1},
    {"relabel",      required_argument,      0,'L'},
    {"lut",          required_argument,      0,'L'},
    {"lut-string",   required_argument,      0,'s'},
    {"aa",           required_argument,      0,'a'},
    {"baa",          no_argument,       &opt.baa_smooth, 1},
    {"transform-bbox",no_argument,      &opt.transform_bbox, 1},
    {"fill",         required_argument, 0, 'f'},
    {0, 0, 0, 0}
    };
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "vqcl:t:o:u:L:l:s:a:U:f:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
      {
      case 0:
          break;
      case 'v':
        cout << "Version: 0.8" << endl;
        return 0;
      case 'a':
        opt.aa_fwhm=atof(optarg);break;
      case 'l':
        like_f=optarg; break;
      case 't':
        xfm_f=optarg; break;
      case 'o':
        order=atoi(optarg);
        order_was_set=true;
        break;
      case 'u':
        opt.uniformize=atof(optarg);break;
      case 'U':
        opt.unistep=atof(optarg);break;
      case 'f':
        opt.fill=atof(optarg);break;
      case 'L':
        map_f=optarg;break;
      case 's':
        map_str=optarg;break;
      case '?':
                /* getopt_long already printed an error message. */
      default:
                show_usage (argv[0]);
                return 1;
        }
    }

  if(opt.normalize && !order_was_set) //perform nearest neighbour interpolation in this case
    order=0;
  
    if ((argc - optind) < 2) {
            show_usage(argv[0]);
            return 1;
    }
  input_f=argv[optind];
  output_f=argv[optind+1];
  
  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }

  if(!map_f.empty())
  {
    if(verbose)
      std::cout<<"Going to relabel input according to:"<<map_f.c_str()<<std::endl;
    std::ifstream in_map(map_f.c_str());
    
    while(!in_map.eof() && in_map.good())
    {
      int l1,l2;
      in_map>>l1>>l2;
      std::map<int,int>::iterator pos=opt.label_map.find(l1);
      if(pos==opt.label_map.end())
        opt.label_map[l1]=l2;
      else if(opt.label_map[l1]!=l2 and verbose)
      {
        std::cerr<<"Warning: label "<<l1<<" mapped more then once!"<<std::endl;
        std::cerr<<"Previous map:"<<(*pos).second<<" New map:"<<l2<<std::endl;
      }
    }
  }
  
  if(!map_str.empty())
  {
    if(verbose)
      std::cout<<"Going to relabel input according to:"<<map_str.c_str()<<std::endl;

    const char* delim1=";";
    char *saveptr1;
    char *_str=strdup(map_str.c_str());
    
    for(char *tok1=strtok_r(_str,delim1,&saveptr1);tok1;tok1=strtok_r(NULL,delim1,&saveptr1))
    {
      int l1,l2;
      if(sscanf(tok1,"%d %d",&l1,&l2)!=2)
      {
        std::cerr<<"Unrecognized lut string:"<<tok1<<std::endl;
      } else {
        std::map<int,int>::iterator pos=opt.label_map.find(l1);
        if(pos==opt.label_map.end())
          opt.label_map[l1]=l2;
        else if(opt.label_map[l1]!=l2 and verbose)
        {
          std::cerr<<"Warning: label "<<l1<<" mapped more then once!"<<std::endl;
          std::cerr<<"Previous map:"<<(*pos).second<<" New map:"<<l2<<std::endl;
        }
      }
    }
    free(_str);
  }

  try
  {
#if ( ITK_VERSION_MAJOR < 4 )  
    itk::RegisterMincIO();
#endif    
    //try to figure out what we have got
    IOBasePointer io = itk::ImageIOFactory::CreateImageIO(input_f.c_str(), itk::ImageIOFactory::ReadMode );
    
    if(!io)
      throw itk::ExceptionObject("Unsupported image file type");
    
    io->SetFileName(input_f.c_str());
    io->ReadImageInformation();

    size_t nd = io->GetNumberOfDimensions();
    size_t nc = io->GetNumberOfComponents();
   
    if(verbose)
    {
      if(opt.uniformize!=0.0)
        std::cout<<"Making uniform sampling, step size="<<opt.uniformize<<std::endl;
      else if(opt.unistep) 
        std::cout<<"Making same step size, new step size="<<opt.unistep<<std::endl;
      else if(opt.normalize)
        std::cout<<"Performing direction cosine normalization"<<std::endl;
    }
    
    if( nc==1 && nd==3 ) //3D image, simple case
    {
      if(labels)
      {
        if(order==0 || !order_was_set)
        {
          //creating the interpolator
          NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
          if(opt.store_byte)
            resample_image<Int3DImage,Byte3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
          else if(opt.store_short)
            resample_image<Int3DImage,Short3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
          else
            resample_image<Int3DImage,Int3DImage,NNInterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
        } else { //using slow algorithm
          InterpolatorType::Pointer interpolator = InterpolatorType::New();
          interpolator->SetSplineOrder(order);
          
          if(opt.store_byte)
            resample_label_image<Int3DImage,Float3DImage,Byte3DImage,InterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
          else if(opt.store_short)
            resample_label_image<Int3DImage,Float3DImage,Short3DImage,InterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
          else
            resample_label_image<Int3DImage,Float3DImage,Int3DImage,InterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
        }
      } else {
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetSplineOrder(order);
        resample_image<Float3DImage,Float3DImage,InterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
      }
     } else if( nd==3 && nc==3 )  { // deal with this case as vector image 
       VectorInterpolatorType::Pointer interpolator = VectorInterpolatorType::New();

       //interpolator->SetSplineOrder(order); // TODO: fix this for ITKv5

       resample_vector_image<Vector3DImage,Vector3DImage,VectorInterpolatorType>(io,output_f,xfm_f,like_f,interpolator,opt);
    } else {
      throw itk::ExceptionObject("This number of dimensions is not supported currently");
    }
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

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; hl c++*/
