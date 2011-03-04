#ifndef _MINC_HELPERS_H_
#define _MINC_HELPERS_H_

#include <complex>
#include <vector>
#include <algorithm>

#include <itkArray.h>
#include <itkOrientedImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>

#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include <minc_io_fixed_vector.h>

namespace minc
{
  
  typedef unsigned char minc_mask_voxel;
  
	typedef float voxel_type;
	const int volume_dimensions = 3;
  
  typedef itk::Vector<float,volume_dimensions>    def_vector;
  typedef itk::Point<double,volume_dimensions>    tag_point;
  typedef std::vector<tag_point>  tag_points;
  
	typedef std::complex < voxel_type > complex;

  
	typedef itk::OrientedImage < complex, volume_dimensions > image3d_complex;
	typedef itk::OrientedImage < voxel_type, volume_dimensions > image3d;
	typedef itk::OrientedImage < minc_mask_voxel, volume_dimensions > mask3d;
	typedef itk::OrientedImage < def_vector, volume_dimensions > def3d;
	
	typedef itk::ImageRegionIteratorWithIndex < image3d > image3d_iterator;
	typedef itk::ImageRegionConstIteratorWithIndex < image3d > image3d_const_iterator;
	
	typedef itk::ImageRegionIteratorWithIndex < mask3d > mask3d_iterator;
	typedef itk::ImageRegionConstIteratorWithIndex < mask3d > mask3d_const_iterator;
	
	typedef itk::ImageRegionIteratorWithIndex < def3d > def3d_iterator;
	typedef itk::ImageRegionConstIteratorWithIndex < def3d > def3d_const_iterator;
	
	typedef itk::ImageRegionIteratorWithIndex < image3d_complex > image3d_complex_iterator;
	typedef itk::ImageRegionConstIteratorWithIndex < image3d_complex > image3d_complex_const_iterator;

  //! find a maximum of elements
  template<class T> float v_max(const T& c)
  {
    float s=std::numeric_limits < float >::min ();;
    for(unsigned int i=0;i<3;i++)
      if(c[i]>s) s=c[i];
    return s;
  }

  //! find a minimum of elements
  template<class T> float v_min(const T& c) 
  {
    float s=std::numeric_limits < float >::max ();;
    for(unsigned int i=0;i<3;i++)
      if(c[i]<s) s=c[i];
    return s;
  }

  //! allocate volume of the same dimension,spacing and origin
  template<class T,class S> void allocate_same(T image,S sample)
  {
		image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
		image->SetBufferedRegion(sample->GetLargestPossibleRegion());
		image->SetRequestedRegion(sample->GetLargestPossibleRegion());
		image->SetSpacing( sample->GetSpacing() );
		image->SetOrigin ( sample->GetOrigin() );
    image->SetDirection(sample->GetDirection());
		image->Allocate();
  }
  
  //! allocate volume of the same dimension,spacing and origin and do nearest neighbour resampling
  template<class T,class S,class E,class D> void nearest_resample_like(T dst,S sample,E src, const D& def)
  {
    allocate_same(dst,sample);    
		itk::ImageRegionIteratorWithIndex<typename T::ObjectType> it(dst, dst->GetLargestPossibleRegion());
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      tag_point p;
      typename E::ObjectType::IndexType idx;
      dst->TransformIndexToPhysicalPoint(it.GetIndex(),p);
      if(src->TransformPhysicalPointToIndex(p,idx))
      {
        it.Value()=src->GetPixel(idx);
      }else{
        it.Value()=def;
      }
    }
  }
  
  //! check if volumes have the same dimensions, spacing and origin
  template<class T,class S> bool check_same(T image,S sample)
  {
		return
      (image->GetLargestPossibleRegion() == sample->GetLargestPossibleRegion()) &&
      (image->GetSpacing() == sample->GetSpacing()) &&
      (image->GetOrigin() == sample->GetOrigin()) &&
      (image->GetDirection().GetVnlMatrix()  == sample->GetDirection().GetVnlMatrix()); // carefull here! , maybe we should calculate some kind of difference here ?
      // this is warkaround a bug in itk
  }
  
  //! allocate volume
  //! \param[out] image - volume to allocate
  //! \param dims - dimensions (voxels)
  //! \param spacing - volume spacing (mm)
  //! \param origin  - volume origin (mm)
  template<class T> void allocate_image3d(T image, 
      const itk::Array<unsigned int> &dims, 
      const itk::Array<double>& spacing, 
      const itk::Array<double>& origin)
	{
		image3d_complex::SizeType  imageSize3D = {{ dims[0], dims[1], dims[2]}};
		image3d_complex::IndexType startIndex3D = { {0, 0, 0}};
		image3d_complex::RegionType region;
		region.SetSize  (imageSize3D);
		region.SetIndex (startIndex3D);
		image->SetLargestPossibleRegion (region);
		image->SetBufferedRegion (region);
		image->SetRequestedRegion (region);
		image->SetSpacing( spacing );
		image->SetOrigin( origin );
		image->Allocate ();
	}
  
  //! allocate volume
  //! \param[out] image - volume to allocate
  //! \param dims - dimensions (voxels)
  //! \param spacing - volume spacing (mm)
  //! \param origin  - volume origin (mm)
  template<class T> void allocate_image3d(T image, 
      const itk::Vector<unsigned int,3> &dims, 
      const itk::Vector<double,3>& spacing, 
      const itk::Vector<double,3>& origin)
	{
		image3d_complex::SizeType  imageSize3D = {{ dims[0], dims[1], dims[2]}};
		image3d_complex::IndexType startIndex3D = { {0, 0, 0}};
		image3d_complex::RegionType region;
		region.SetSize  (imageSize3D);
		region.SetIndex (startIndex3D);
		image->SetLargestPossibleRegion (region);
		image->SetBufferedRegion (region);
		image->SetRequestedRegion (region);
		image->SetSpacing( spacing );
		image->SetOrigin( origin.GetDataPointer () );
		image->Allocate ();
	}
  
  //! allocate volume
  //! \param[out] image - volume to allocate
  //! \param dims - dimensions (voxels)
  //! \param spacing - volume spacing (mm)
  //! \param origin  - volume origin (mm)
  template<class T> void allocate_image3d(T image, 
      const fixed_vec<3, unsigned int>&dims, 
      const fixed_vec<3, double>& spacing=fixed_vec<3, double>(1.0) , 
      const fixed_vec<3, double>& origin=fixed_vec<3, double>(0.0))
	{
		image3d_complex::SizeType  imageSize3D = {{ dims[0], dims[1], dims[2]}};
		image3d_complex::IndexType startIndex3D = { {0, 0, 0}};
		image3d_complex::RegionType region;
		region.SetSize  (imageSize3D);
		region.SetIndex (startIndex3D);
		image->SetLargestPossibleRegion (region);
		image->SetBufferedRegion (region);
		image->SetRequestedRegion (region);
		image->SetSpacing( spacing.c_buf() );
		image->SetOrigin( origin.c_buf() );
		image->Allocate ();
	}

  
  inline image3d::SizeType operator/= (image3d::SizeType & s, int d)
  {
    s[0] /= d;
    s[1] /= d;
    s[2] /= d;
    return s;
  }

  inline image3d::SizeType operator*= (image3d::SizeType & s, int d)
  {
    s[0] *= d;
    s[1] *= d;
    s[2] *= d;
    return s;
  }
  
// a helper function for minc reading
  template <class T> typename T::Pointer load_minc(const char *file)
  {
    typedef itk::MincImageIO ImageIOType;
    ImageIOType::Pointer minc2ImageIO = ImageIOType::New();
     
    typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
    
    reader->SetFileName(file);
    reader->SetImageIO( minc2ImageIO );
    reader->Update();
    
    return reader->GetOutput();
  }

  
  void set_minc_storage_type(itk::Object* image,nc_type datatype,bool is_signed);
  void copy_metadata(itk::Object* dst,itk::Object* src);
  void append_history(itk::Object* dst,const std::string& history);
   
  // a helper function for minc writing
  template <class T> void save_minc(const char *file,const T& img)
  {
    typedef itk::MincImageIO ImageIOType;
    ImageIOType::Pointer minc2ImageIO = ImageIOType::New();
     
    typename itk::ImageFileWriter< T >::Pointer writer = itk::ImageFileWriter<T>::New();
    writer->SetFileName(file);
    writer->SetImageIO( minc2ImageIO );
    writer->SetInput( img );
    writer->Update();
  } 
  
  // a helper function for minc writing
  template <class T> void save_minc(const char *file,const T* img)
  {
    typedef itk::MincImageIO ImageIOType;
    ImageIOType::Pointer minc2ImageIO = ImageIOType::New();
     
    typename itk::ImageFileWriter< T >::Pointer writer = itk::ImageFileWriter<T>::New();
    writer->SetFileName(file);
    writer->SetImageIO( minc2ImageIO );
    writer->SetInput( img );
    writer->Update();
  } 
  
  //! allocate volume with spacing 1mm and origin 0,0,0
  //! \param[out] image - volume to allocate
  //! \param dims - dimensions (voxels)
  template<class T> void allocate_image3d(T image, const itk::Size<3> &dims)
  {
    allocate_image3d(image,IDX<unsigned int>(dims[0],dims[1],dims[2]));
  }

  //! calculate volume min and max
  int get_image_limits(image3d::Pointer, voxel_type &min, voxel_type &max);
  //! calculate volume min and max
  int get_image_limits(def3d::Pointer, voxel_type &min, voxel_type &max);
  //! calculate volume min and max
  int get_image_limits(mask3d::Pointer img, voxel_type &min,voxel_type &max);
  
  //! store tags to minc tags file
  void write_tags(const tag_points& tags, const char * file);
  //! store two sets of tags to minc tags file
  void write_2tags(const tag_points& tags,const tag_points& tags2, const char * file);
  
  //! store tags and labels
  void write_tags(const tag_points& tags,const std::vector<int>& labels, const char * file);
  //! store tags and values
  void write_tags(const tag_points& tags,const std::vector<float>& values, const char * file);
  //! store tags and values
  void write_tags(const tag_points& tags,const std::vector<double>& values, const char * file);
  
  
  //! read tags from the minc tag file
  void read_tags(tag_points& tags, const char * file,int vol=1);
  
  //! read tags and labels from minc tag file
  void read_tags(tag_points& tags, std::vector<int>& labels, const char * file, int vol=1);
  
  //! read tags and labels from minc tag file
  void read_tags(tag_points& tags, std::vector<float>& labels, const char * file, int vol=1);
  
  //! read tags and labels from minc tag file
  void read_tags(tag_points& tags, std::vector<double>& labels, const char * file, int vol=1);

  //! read tags and labels from minc tag file
  void read_tags(tag_points& tag1, tag_points& tag2,std::vector<double>& labels, const char * file);
  
  
  //! read array from the text file
  void load_parameters(const char *file,itk::Array<double> &param);
  //! save array to the text file
  void save_parameters(const char *file,const itk::Array<double> &param);
  //! read array from the text file, up to components
  void load_parameters(const char *file,itk::Array<double> &param,int no);
  //! make sure that all voxels are 0 or 1
  void normalize_mask(mask3d::Pointer img);
  
  void read_linear_xfm(const char *xfm,itk::Matrix<double,3,3>& rot, itk::Vector<double,3>& tran );
  void read_linear_xfm(const char *xfm,itk::Matrix<double,2,2>& rot, itk::Vector<double,2>& tran );
  
  void write_linear_xfm(const char *xfm,const itk::Matrix<double,3,3>& rot,const itk::Vector<double,3>& tran);
  void write_linear_xfm(const char *xfm,const itk::Matrix<double,2,2>& rot,const itk::Vector<double,2>& tran);
  
  void write_combined_xfm(const char *xfm,const itk::Matrix<double,3,3>& rot,const itk::Vector<double,3>& tran, const char *grid );
  void write_combined_xfm(const char *xfm,const itk::Matrix<double,2,2>& rot,const itk::Vector<double,2>& tran, const char *grid );
  
  void write_combined_xfm(const char *xfm,const char *grid, const itk::Matrix<double,3,3>& rot,const itk::Vector<double,3>& tran );
  void write_combined_xfm(const char *xfm,const char *grid, const itk::Matrix<double,2,2>& rot,const itk::Vector<double,2>& tran );
  
  void write_nonlinear_xfm(const char *xfm,const char *grid);
};

#endif //_MINC_HELPERS_H_
