#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkImageIOBase.h>
#include <itkFlipImageFilter.h>
#include <itkMetaDataObject.h>
#include <itkMetaDataDictionary.h>
#include <itkCastImageFilter.h>
#include <itkDiffusionTensor3D.h>

#if ( ITK_VERSION_MAJOR > 3 ) 
  #include "itk4MincHelpers.h"
#else
  #include "itkMincImageIOFactory.h"
  #include "itkMincImageIO.h"
  #include <time_stamp.h>    // for creating minc style history entry
  #include "itkMincHelpers.h"
#endif

#include <getopt.h>
#include <unistd.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>


typedef std::vector<double> double_vector;

//TODO: add support for 6-component symmetric matrix file type

void show_usage (const char *name)
{
  std::cerr  
      << "Usage: "<<name<<" <input> <output> " << std::endl
      << "--clobber clobber output files"<<std::endl
      << "--verbose be verbose"<<std::endl
      << "--show-meta show meta information (if present)"<<std::endl
      << " Spatial transformation flags:"<<std::endl
      << "--inv-x invert X axis"<<std::endl
      << "--inv-y invert Y axis"<<std::endl
      << "--inv-z invert Z axis"<<std::endl
      << "--center set origin to the center of the image"<<std::endl
      << " DWI related flags "<<std::endl
      << "--dwi - assume that we are dealing with DWI scan"<<std::endl
      << "--dti - assume that we are dealing with DTI information (or DWI scan)"<<std::endl
      << "--use-b-matrix - convert b-matrix (if present) into DWI gradient directions and b-value, implies DWI"<<std::endl
      << "--minc-to-nrrd Convert minc style to nrrd style "<<std::endl
      << "--nrrd-to-minc Convert nrrd style to minc style "<<std::endl
      << " Data conversion flags: "<<std::endl
      << "--char   - cast data as signed char"<<std::endl
      << "--byte   - cast data as unsigned char"<<std::endl
      << "--short  - cast data as signed short"<<std::endl
      << "--ushort - cast data as unsigned short"<<std::endl
      << "--int    - cast data as signed int"<<std::endl
      << "--uint   - cast data as unsigned int"<<std::endl
      << "--float  - cast data as float"<<std::endl
      << "--double - cast data as double"<<std::endl
      << " MINC file format related flags:" << std::endl
      << "--mbyte   - save pixel value using signed char"<<std::endl
      << "--mshort  - save pixel value using short"<<std::endl
      << "--mint    - save pixel value using int"<<std::endl
      << "--mfloat  - save pixel value using float"<<std::endl
      << "--mdouble - save pixel value using double"<<std::endl
      
;
}

typedef itk::ImageIOBase IOBase;
typedef itk::SmartPointer<IOBase> IOBasePointer;

class ImageConverterBase
{
protected:
  bool assume_dti;
  bool use_b_matrix;
  bool dwi_flip_z;
  std::string history;
  bool verbose;
  bool inv_x;
  bool inv_y;
  bool inv_z;
  bool center;
  bool show_meta;
  std::string  minc_type;
  bool minc_to_nrrd;
  bool nrrd_to_minc;
public:
  void setup(bool _verbose=false,
              bool _assume_dti=false,
              bool _use_b_matrix=false,
              bool _dwi_flip_z=false,
              bool _inv_x=false,bool _inv_y=false,bool _inv_z=false, 
              bool _center=false, 
              bool _show_meta=false,
              const std::string& _history="",
              const std::string& _minc_type="",
              bool _minc_to_nrrd=false,
              bool _nrrd_to_minc=false
            )
  {
    assume_dti=_assume_dti;
    use_b_matrix=_use_b_matrix;
    dwi_flip_z=_dwi_flip_z;
    history=_history;
    verbose=_verbose;
    inv_x=_inv_x;
    inv_y=_inv_y;
    inv_z=_inv_z;
    center=_center;
    show_meta=_show_meta;
    minc_type=_minc_type;
    minc_to_nrrd=_minc_to_nrrd;
    nrrd_to_minc=_nrrd_to_minc;
  }
  
  virtual void load_and_save_image(IOBase* base,const char *fname,itk::ImageIOBase::IOComponentType oct)=0;

  virtual ~ImageConverterBase()
  {}

  void convert_meta_minc_to_nrrd(itk::MetaDataDictionary& dict)
  {
    // let's try converting DTI-related meta-information
    itk::Array<double> bvalues,direction_x,direction_y,direction_z;
    
    if( itk::ExposeMetaData< itk::Array<double> >( dict , "acquisition:bvalues",bvalues) &&
        itk::ExposeMetaData< itk::Array<double> >( dict , "acquisition:direction_x",direction_x) &&
        itk::ExposeMetaData< itk::Array<double> >( dict , "acquisition:direction_y",direction_y) &&
        itk::ExposeMetaData< itk::Array<double> >( dict , "acquisition:direction_z",direction_z))
    {
      //We have got MINC-style DTI metadata
      if(bvalues.size()!=direction_x.size() || 
         bvalues.size()!=direction_y.size() ||
         bvalues.size()!=direction_z.size() )
      {
        std::cerr<<"WARNING: Different number of components of DTI directions"<<std::endl
                <<" Skipping DTI metadata conversion!"<<std::endl;
      } else {
        double max_b_value=*std::max_element(bvalues.begin(),bvalues.end());
        if(verbose)
          std::cout<<"Found maximum B-value:"<<max_b_value<<std::endl;
        
        itk::EncapsulateMetaData<std::string>( dict,
          "modality",
          "DWMRI");      
        
        //TODO: find out if we can get this info from minc 
        // using identity for now
        std::vector< std::vector< double > > measurement_frame(3,std::vector< double >(3));
        measurement_frame[0][0]=1.0;
        measurement_frame[1][1]=1.0;
        measurement_frame[2][2]=1.0;
        
        
        itk::EncapsulateMetaData< std::vector< std::vector< double > > >( dict,
          "NRRD_measurement frame",
          measurement_frame);      

        
        //TODO: deal with NRRD_thicknesses ?
        
        std::ostringstream ossKey;
        ossKey<<std::setw(9) << std::setiosflags(std::ios::fixed)
          << std::setprecision(6) << std::setiosflags(std::ios::right) << max_b_value;
        
          
        itk::EncapsulateMetaData<std::string>( dict,
          "NRRD_centerings[0]","cell");
        itk::EncapsulateMetaData<std::string>( dict,
          "NRRD_centerings[1]","cell");
        itk::EncapsulateMetaData<std::string>( dict,
          "NRRD_centerings[2]","cell");
        
        itk::EncapsulateMetaData<std::string>( dict,
          "NRRD_kinds[0]","space");
        itk::EncapsulateMetaData<std::string>( dict,
          "NRRD_kinds[1]","space");
        itk::EncapsulateMetaData<std::string>( dict,
          "NRRD_kinds[2]","space");

        
        itk::EncapsulateMetaData<std::string>( dict,
          "DWMRI_b-value",
          ossKey.str());
        
        
        size_t zero_b_value_cnt=0;
        for(size_t i=0;i<bvalues.size();i++)
        {
          double direction[3]={direction_x[i],direction_y[i],direction_z[i]};
          double bval=bvalues[i];
          double b_scale=1.0;
          
          if(fabs(bval)<1e-3) 
            zero_b_value_cnt++;
          else
          {
            b_scale=::sqrt(bval/max_b_value);
            
            direction[0]*=b_scale;
            direction[1]*=b_scale;
            direction[2]*=(dwi_flip_z&&!use_b_matrix)?-b_scale:b_scale;
          }
          
          std::ostringstream ossKey2;
          ossKey2 << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << i;

          std::ostringstream ossMetaString;
          ossMetaString << std::setw(9) << std::setiosflags(std::ios::fixed)
            << std::setprecision(6) << std::setiosflags(std::ios::right)
            << direction[0]
          << "    "
            << std::setw(9) << std::setiosflags(std::ios::fixed)
            << std::setprecision(6) << std::setiosflags(std::ios::right)
            << direction[1]
          << "    "
            << std::setw(9) << std::setiosflags(std::ios::fixed)
            << std::setprecision(6) << std::setiosflags(std::ios::right)
            << direction[2];

          // std::cout<<ossKey.str()<<ossMetaString.str()<<std::endl;
          itk::EncapsulateMetaData<std::string>( dict,
            ossKey2.str(), ossMetaString.str() );
        }
        if(verbose)
          std::cout<<"Found "<<zero_b_value_cnt<<" Zero b-value components"<<std::endl;
      }
    }
  }

  void convert_meta_nrrd_to_minc(itk::MetaDataDictionary& dict)
  {
    //1. get b-value
    std::string b_value_,gradient_value_;
    double bval=0.0,vect3d[3];
    double_vector bvalues,direction_x,direction_y,direction_z;
    
    //remove B-matrix, it is of no use now
    //if(dict.HasKey( "acquisition:b_matrix" ))
      
    
    if(itk::ExposeMetaData<std::string>( dict,"DWMRI_b-value", b_value_))
    {
      int i=0;
      while(true)
      {
          bval=atof(b_value_.c_str());
          
          std::ostringstream ossKey;
          ossKey << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << i++;
          if(itk::ExposeMetaData<std::string>( dict,ossKey.str(), gradient_value_))
          {
            std::istringstream iss(gradient_value_);
            iss >> vect3d[0] >> vect3d[1] >> vect3d[2];
            
            //? normalize to unit vector and modulate bvalue
            double b_scale=vect3d[0]*vect3d[0]+vect3d[1]*vect3d[1]+vect3d[2]*vect3d[2];
            double vmag=::sqrt(b_scale);
            
            if(vmag<0.1) //Zero
            {
              bvalues.push_back(0.0); 
              direction_x.push_back(0.0);
              direction_y.push_back(0.0);
              direction_z.push_back(0.0);
            } else {
              bvalues.push_back(bval*b_scale);
              direction_x.push_back(vect3d[0]/vmag);
              direction_y.push_back(vect3d[1]/vmag);
              direction_z.push_back(((dwi_flip_z&&!use_b_matrix)?-vect3d[2]:vect3d[2])/vmag);
            }
          } else {
            break;
          }
      }
      
      if(!bvalues.empty())
      {
        itk::EncapsulateMetaData<double_vector>( dict , "acquisition:bvalues",bvalues);
        itk::EncapsulateMetaData<double_vector>( dict , "acquisition:direction_x",direction_x);
        itk::EncapsulateMetaData<double_vector>( dict , "acquisition:direction_y",direction_y);
        itk::EncapsulateMetaData<double_vector>( dict , "acquisition:direction_z",direction_z);
        itk::EncapsulateMetaData<double>( dict ,        "acquisition:b_value", bval);

        // cleanup NRRD style meta information
        // unfortunately it is impossible to delete entries from metadata :(
        itk::EncapsulateMetaData<std::string>( dict , "DWMRI_b-value","");
        for(int i=0;i<bvalues.size();i++)
        {
          std::ostringstream ossKey;
          ossKey << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << i;
          itk::EncapsulateMetaData<std::string>( dict ,ossKey.str(),"");
        }
      }
    }
    
    
#if ( ITK_VERSION_MAJOR < 4 ) 
    // force time to be fastes varying dimension to simplify using xdisp
    std::vector<std::string> dimorder;

    dimorder.push_back(MItime); //fastest varying
    
    dimorder.push_back(MIzspace);
    dimorder.push_back(MIyspace);
    dimorder.push_back(MIxspace);
    
    itk::EncapsulateMetaData(dict,"dimorder", dimorder);
#endif    

  }



  double decompose_b_matrix(const double_vector& b_matrix, 
                          double_vector & bvalues,
                          double_vector & direction_x, 
                          double_vector & direction_y,
                          double_vector & direction_z)
  {
    size_t elements=b_matrix.size();
    if(elements%6) //should be divisible by 6 
      std::cerr<<"Number of elements in b_matrix is not divisible by 6 , probably an error!"<<std::endl;
    elements/=6;
    
    bvalues.resize(elements);
    direction_x.resize(elements);
    direction_y.resize(elements);
    direction_z.resize(elements);

    vnl_matrix_fixed<double, 3, 3> bMatrix;
    double max_bval=0;
    
    for(size_t i=0;i<elements;i++)
    {
      //Code from DicomToNrrdConverter
      // UNC comments: The principal eigenvector of the bmatrix is to be extracted as
      // it's the gradient direction and trace of the matrix is the b-value

      // UNC comments: Fill out the 3x3 bmatrix with the 6 components read from the 
      // DICOM header.
      bMatrix[0][0] = b_matrix[i*6+0];
      bMatrix[0][1] = b_matrix[i*6+1];
      bMatrix[0][2] = b_matrix[i*6+2];
      bMatrix[1][1] = b_matrix[i*6+3];
      bMatrix[1][2] = b_matrix[i*6+4];
      bMatrix[2][2] = b_matrix[i*6+5];
      bMatrix[1][0] = bMatrix[0][1];
      bMatrix[2][0] = bMatrix[0][2];
      bMatrix[2][1] = bMatrix[1][2];

      double bvalue = bMatrix[0][0] + bMatrix[1][1] + bMatrix[2][2];
      
      if(bvalue>max_bval) max_bval=bvalue;
      
      // UNC comments: The b-value si the trace of the bmatrix
      // UNC comments: Even if the bmatrix is null, the svd decomposition set the 1st eigenvector
      // to (1,0,0). So we force the gradient direction to 0 if the bvalue is null
      if(bvalue == 0.0 ) //TODO: threshold?
      {
        direction_x[i] = 0.0;
        direction_y[i] = 0.0;
        direction_z[i] = 0.0;
        bvalues[i] = 0.0;
      } else {
        // UNC comments: Computing the decomposition
        vnl_svd<double> svd(bMatrix.as_ref());

        // UNC comments: Extracting the principal eigenvector i.e. the gradient direction
        
        direction_x[i] = svd.U(0,0);
        direction_y[i] = svd.U(1,0);
        direction_z[i] = dwi_flip_z?svd.U(2,0):svd.U(2,0);
        bvalues[i] = bvalue;
      }
    }
    return max_bval;
  }
  
  void print_metadata(itk::MetaDataDictionary &thisDic)
  {
    //let's write some meta information if there is any 
    for(itk::MetaDataDictionary::ConstIterator it=thisDic.Begin();it!=thisDic.End();++it)
    {
      itk::MetaDataObjectBase *bs=(*it).second;
      itk::MetaDataObject<std::string> * str=dynamic_cast<itk::MetaDataObject<std::string> *>(bs);
      if(str)
        std::cout<<(*it).first.c_str()<<" = "<< str->GetMetaDataObjectValue().c_str()<<std::endl;
      else
        std::cout<<(*it).first.c_str()<<" type: "<< typeid(*bs).name()<<std::endl;
    }
  }
  
  void convert_bmatrix(itk::MetaDataDictionary &thisDic)
  {
    double_vector bmatrix;
    
    if( itk::ExposeMetaData<double_vector>( thisDic , "acquisition:b_matrix",bmatrix))
    {
      double_vector bvalues,direction_x,direction_y,direction_z;
      double bval=decompose_b_matrix(bmatrix,bvalues,direction_x,direction_y,direction_z);
      
      if(!bvalues.empty())
      {
        itk::EncapsulateMetaData<double_vector>( thisDic , "acquisition:bvalues"    ,bvalues);
        itk::EncapsulateMetaData<double_vector>( thisDic , "acquisition:direction_x",direction_x);
        itk::EncapsulateMetaData<double_vector>( thisDic , "acquisition:direction_y",direction_y);
        itk::EncapsulateMetaData<double_vector>( thisDic , "acquisition:direction_z",direction_z);
        itk::EncapsulateMetaData<double>(        thisDic , "acquisition:b_value"    ,bval);
      }
    }
  }
  
  void convert_metadata(itk::MetaDataDictionary &thisDic)
  {
    if(show_meta)
     print_metadata(thisDic);
    
    if( thisDic.HasKey( "acquisition:b_matrix" ) && use_b_matrix )
    { 
      convert_bmatrix(thisDic);
    }
    
    //making sure that all vectors contain the same number of parameters (just in case)
    if( thisDic.HasKey( "acquisition:bvalues"    ) &&
        thisDic.HasKey( "acquisition:direction_x") &&
        thisDic.HasKey( "acquisition:direction_y") &&
        thisDic.HasKey( "acquisition:direction_z") && 
        !nrrd_to_minc)
    {
      if(verbose)
        std::cout<<"Converting MINC-style DWI to NRRD style"<<std::endl;
      convert_meta_minc_to_nrrd(thisDic);
    } else if ( thisDic.HasKey("DWMRI_b-value")  && !minc_to_nrrd  ) {
      //We have got NRRD-style DTI metadata
      if(verbose)
        std::cout<<"Converting NRRD-style DWI to MINC style"<<std::endl;
      convert_meta_nrrd_to_minc(thisDic);
    } else if( nrrd_to_minc || minc_to_nrrd) {
      std::cerr<<"ERROR: requested DWI headers are missing!"<<std::endl;
    }
  }
};

template<class TInputImage> class ImageConverter: public ImageConverterBase
{
protected:
  
  typedef TInputImage InputImageType;
  typedef typename InputImageType::PixelType InputImagePixelType;

public:

  ~ImageConverter()
  {

  }

  ImageConverter()
  {
    setup();
  }

  template<class TOutputImage> void _load_and_save_image(IOBase* base,const char *fname)
  {
    typename itk::ImageFileReader<TInputImage >::Pointer reader = itk::ImageFileReader<TInputImage >::New();
    typename itk::FlipImageFilter<TInputImage >::Pointer flip=itk::FlipImageFilter<TInputImage >::New();
    
    
    reader->SetImageIO(base);
    reader->SetFileName(base->GetFileName());
    reader->Update();
    
    typename TInputImage::Pointer img=reader->GetOutput();
    itk::MetaDataDictionary &thisDic=img->GetMetaDataDictionary();

    convert_metadata(thisDic);
    
    // let's analyze direction cosines, and convert image into something resambling neurological notation
    if(verbose)
    {
      std::cout<<"Dimensions:"<<img->GetLargestPossibleRegion().GetSize()<<" ";
      std::cout<<"Origin:"<<img->GetOrigin()<<std::endl;
      std::cout<<"Steps:"<<img->GetSpacing()<<std::endl;
      std::cout<<"Directions:["<<img->GetDirection()<<"]"<<std::endl;
    }
    
    if(inv_x||inv_y||inv_z)
    {
      typename itk::FlipImageFilter<TInputImage >::FlipAxesArrayType arr;
      arr[0]=inv_x;
      arr[1]=inv_y;
      arr[2]=inv_z;
      flip->SetFlipAxes(arr);
      flip->SetInput(img);
      flip->Update();
      img=flip->GetOutput();
    }
    
    if(center)//move origin to the center of the image
    {
      typename TInputImage::RegionType r=img->GetLargestPossibleRegion();
      //std::vector<double> corner[3];
      
      typename TInputImage::IndexType idx;
      typename TInputImage::PointType c;
      
      idx[0]=r.GetIndex()[0]+r.GetSize()[0]/2.0;
      idx[1]=r.GetIndex()[1]+r.GetSize()[1]/2.0;
      idx[2]=r.GetIndex()[2]+r.GetSize()[2]/2.0;
      
      img->TransformIndexToPhysicalPoint(idx,c);
      
      typename TInputImage::PointType org=img->GetOrigin();
      
      org[0]-=c[0];
      org[1]-=c[1];
      org[2]-=c[2];
      
      img->SetOrigin(org);
    }
    
    if(!history.empty())
      minc::append_history(img,history);
      
    typename itk::CastImageFilter< TInputImage, TOutputImage >::Pointer cast=itk::CastImageFilter< TInputImage, TOutputImage >::New();
    
    cast->SetInput(img);
    
    if(verbose)
      std::cout<<"Writing "<<fname<<"..."<<std::endl;
    
    typename itk::ImageFileWriter< TOutputImage >::Pointer writer = itk::ImageFileWriter<TOutputImage >::New();
    writer->SetFileName(fname);
    cast->Update();
    
    minc::copy_metadata(cast->GetOutput(),img);
    
    if(!minc_type.empty()) 
      minc::set_minc_storage_type(cast->GetOutput(),minc_type);
    
    if(verbose)
    {
      std::cout<<"minc_type="<<minc_type<<std::endl;
    }
    
    writer->SetInput( cast->GetOutput() );
    writer->Update();
  }

  virtual void load_and_save_image(IOBase* io,const char *fname,itk::ImageIOBase::IOComponentType oct)
  {
    switch(oct)
    {
      case itk::ImageIOBase::UCHAR:
        this->_load_and_save_image<itk::VectorImage<unsigned char, InputImageType::ImageDimension> >(io,fname);
        break;
      case itk::ImageIOBase::CHAR:
        this->_load_and_save_image<itk::VectorImage<char, InputImageType::ImageDimension> >(io,fname);
        break;
      case itk::ImageIOBase::USHORT:
        _load_and_save_image<itk::VectorImage<unsigned short, InputImageType::ImageDimension> >(io,fname);
        break;
      case itk::ImageIOBase::SHORT:
        _load_and_save_image<itk::VectorImage<short, InputImageType::ImageDimension> >(io,fname);
        break;
      case itk::ImageIOBase::INT:
        _load_and_save_image<itk::VectorImage<int, InputImageType::ImageDimension> >(io,fname);
        break; 
      case itk::ImageIOBase::UINT:
        _load_and_save_image<itk::VectorImage<unsigned int, InputImageType::ImageDimension> >(io,fname);
        break; 
      case itk::ImageIOBase::FLOAT:
        _load_and_save_image<itk::VectorImage<float, InputImageType::ImageDimension> >(io,fname);
        break; 
      case itk::ImageIOBase::DOUBLE:
        _load_and_save_image<itk::VectorImage<double, InputImageType::ImageDimension> >(io,fname);
        break; 
      default:
        itk::ExceptionObject("Unsupported component type");
    }
  }
};

template<class TInputImage> class TensorImageConverter: public ImageConverterBase
{
protected:
  
  typedef TInputImage InputImageType;
  typedef typename InputImageType::PixelType InputImagePixelType;

public:

  ~TensorImageConverter()
  {

  }

  TensorImageConverter()
  {
    setup();
  }

  template<class TOutputImage> void _load_and_save_image(IOBase* base,const char *fname)
  {
    typename itk::ImageFileReader<TInputImage >::Pointer reader = itk::ImageFileReader<TInputImage >::New();
    typename itk::FlipImageFilter<TInputImage >::Pointer flip=itk::FlipImageFilter<TInputImage >::New();
    
    
    reader->SetImageIO(base);
    reader->SetFileName(base->GetFileName());
    reader->Update();
    
    typename TInputImage::Pointer img=reader->GetOutput();
    itk::MetaDataDictionary &thisDic=img->GetMetaDataDictionary();

    convert_metadata(thisDic);
    
    // let's analyze direction cosines, and convert image into something resambling neurological notation
    if(verbose)
    {
      std::cout<<"Dimensions:"<<img->GetLargestPossibleRegion().GetSize()<<" ";
      std::cout<<"Origin:"<<img->GetOrigin()<<std::endl;
      std::cout<<"Steps:"<<img->GetSpacing()<<std::endl;
      std::cout<<"Directions:["<<img->GetDirection()<<"]"<<std::endl;
    }
    
    if(inv_x||inv_y||inv_z)
    {
      typename itk::FlipImageFilter<TInputImage >::FlipAxesArrayType arr;
      arr[0]=inv_x;
      arr[1]=inv_y;
      arr[2]=inv_z;
      flip->SetFlipAxes(arr);
      flip->SetInput(img);
      flip->Update();
      img=flip->GetOutput();
    }
    
    if(center)//move origin to the center of the image
    {
      typename TInputImage::RegionType r=img->GetLargestPossibleRegion();
      //std::vector<double> corner[3];
      
      typename TInputImage::IndexType idx;
      typename TInputImage::PointType c;
      
      idx[0]=r.GetIndex()[0]+r.GetSize()[0]/2.0;
      idx[1]=r.GetIndex()[1]+r.GetSize()[1]/2.0;
      idx[2]=r.GetIndex()[2]+r.GetSize()[2]/2.0;
      
      img->TransformIndexToPhysicalPoint(idx,c);
      
      typename TInputImage::PointType org=img->GetOrigin();
      
      org[0]-=c[0];
      org[1]-=c[1];
      org[2]-=c[2];
      
      img->SetOrigin(org);
    }
    
    if(!history.empty())
        minc::append_history(img,history);
    
    typename itk::CastImageFilter< TInputImage, TOutputImage >::Pointer cast=itk::CastImageFilter< TInputImage, TOutputImage >::New();
    
    cast->SetInput(img);
    
    if(verbose)
      std::cout<<"Writing "<<fname<<"..."<<std::endl;
    
    typename itk::ImageFileWriter< TOutputImage >::Pointer writer = itk::ImageFileWriter<TOutputImage >::New();
    writer->SetFileName(fname);
    cast->Update();
    
    minc::copy_metadata(cast->GetOutput(),img);
    
    if(!minc_type.empty()) 
      minc::set_minc_storage_type(cast->GetOutput(),minc_type);
    
    if(verbose)
    {
      std::cout<<"minc_type="<<minc_type<<std::endl;
    }

    
    writer->SetInput( cast->GetOutput() );
    writer->Update();
  }

  virtual void load_and_save_image(IOBase* io, const char *fname, itk::ImageIOBase::IOComponentType oct)
  {
    switch(oct)
    {
      //default:
      case itk::ImageIOBase::FLOAT:
        _load_and_save_image<itk::Image<itk::DiffusionTensor3D<float>,3 > >(io,fname);
        break; 
      case itk::ImageIOBase::DOUBLE:
        _load_and_save_image<itk::Image<itk::DiffusionTensor3D<double>,3 > >(io,fname);
        break; 
    }
  }
};


int main(int argc,char **argv)
{
  int verbose=0;
  int clobber=0;
  int show_meta=0;
  int inv_x=0,inv_y=0,inv_z=0,center=0;
  int assume_dti=0;
#if ( ITK_VERSION_MAJOR < 4 ) 
  char *_history = time_stamp(argc, argv);
  std::string history=_history;
  free(_history);
#else
  std::string history= minc_timestamp(argc,argv);
#endif
  int use_b_matrix=0;
  int dwi_flip_z=0;
  int nrrd_to_minc=0;
  int minc_to_nrrd=0;
  
  int store_char=0,store_uchar=0,store_short=0,store_ushort=0,store_float=0,store_int=0,store_uint=0,store_double=0;
  int minc_float=0,minc_double=0,minc_byte=0,minc_short=0,minc_int=0;
  std::string minc_type;
  
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"inv-x", no_argument,  &inv_x, 1},
    {"inv-y", no_argument,  &inv_y, 1},
    {"inv-z", no_argument,  &inv_z, 1},
    {"center", no_argument, &center, 1},
    {"show-meta", no_argument, &show_meta, 1},
    {"dti", no_argument, &assume_dti, 1},
    {"dwi", no_argument, &assume_dti, 1},
    {"use-b-matrix",no_argument, &use_b_matrix, 1},
    {"nrrd-to-minc",no_argument, &nrrd_to_minc, 1},
    {"minc-to-nrrd",no_argument, &minc_to_nrrd, 1},
    
    {"float", no_argument,  &store_float,  1},
    {"double", no_argument, &store_double, 1},
    {"byte", no_argument,   &store_uchar,  1},
    {"char", no_argument,   &store_char,   1},
    {"short", no_argument,  &store_short,  1},
    {"ushort", no_argument, &store_ushort, 1},
    {"int", no_argument,    &store_int,    1},
    {"uint", no_argument,   &store_uint,   1},
    
    {"mfloat",  no_argument, &minc_float, 1},
    {"mdouble", no_argument, &minc_double, 1},
    {"mbyte",   no_argument, &minc_byte, 1},
    {"mshort",  no_argument, &minc_short, 1},
    {"mint",    no_argument, &minc_int,   1},
    
    {0, 0, 0, 0}
  };
    
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
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
  if(nrrd_to_minc || minc_to_nrrd)
    dwi_flip_z=1;
  
  try
  {
#if ( ITK_VERSION_MAJOR <4 ) 
    itk::RegisterMincIO();
#endif
    
    if(use_b_matrix) assume_dti=1;
    
    /* READING */
    if( verbose ) 
      std::cout<<"Reading "<<input.c_str()<<"..."<<std::endl;
    
    //try to figure out what we have got
    IOBasePointer io = itk::ImageIOFactory::CreateImageIO(input.c_str(), itk::ImageIOFactory::ReadMode );
    
    if(!io)
      throw itk::ExceptionObject("Unsupported image file type");
    
    io->SetFileName(input.c_str());
    io->ReadImageInformation();

    size_t nd = io->GetNumberOfDimensions();
    size_t nc = io->GetNumberOfComponents();
    itk::ImageIOBase::IOComponentType  ct = io->GetComponentType();
    itk::ImageIOBase::IOComponentType  oct = ct;
    itk::ImageIOBase::IOPixelType      pt = io->GetPixelType();
    
  
    if(store_char)
      oct=itk::ImageIOBase::CHAR;
    else if(store_uchar)
      oct=itk::ImageIOBase::UCHAR;
    else if(store_short)
      oct=itk::ImageIOBase::SHORT;
    else if(store_ushort)
      oct=itk::ImageIOBase::USHORT;
    else if(store_int)
      oct=itk::ImageIOBase::INT;
    else if(store_uint)
      oct=itk::ImageIOBase::UINT;
    else if(store_float)
      oct=itk::ImageIOBase::FLOAT;
    else if(store_double)
      oct=itk::ImageIOBase::DOUBLE;
    
    std::string ct_s = io->GetComponentTypeAsString(ct);
    std::string oct_s = io->GetComponentTypeAsString(oct);
    
    std::string pt_s = io->GetPixelTypeAsString(pt);
    
    if(verbose)
    {
      std::cout<<"dimensions:"<< nd << std::endl
               <<"components:"<< nc << std::endl
               <<"pixel type:"<< pt_s << std::endl
               <<"input type:"<< ct_s.c_str() <<std::endl
               <<"output type:"<<oct_s.c_str()<<std::endl;
    }

    ImageConverterBase *converter=NULL;

    if(nd==3 && nc==6 && assume_dti) //we are dealing with tesor image
    {
      if(verbose) std::cout<<"Writing 3D tensor image..."<<std::endl;
      switch(ct) //TODO: maybe tensors should be stored as float or double only?
      { 
        case itk::ImageIOBase::DOUBLE:
          converter=new TensorImageConverter<itk::Image<itk::DiffusionTensor3D<double>,3> >();
          break; 
        default  : //itk::ImageIOBase::FLOAT 
          converter=new TensorImageConverter<itk::Image<itk::DiffusionTensor3D<float>,3> >();
          break; 
      }
    } else if(nd==3 || assume_dti)
    {
      if(verbose) std::cout<<"Writing 3D image..."<<std::endl;
      switch(ct)
      { 
        case itk::ImageIOBase::UCHAR :
          converter=new ImageConverter<itk::VectorImage<unsigned char, 3> >();
          break;
        case itk::ImageIOBase::CHAR :
          converter=new ImageConverter<itk::VectorImage<char, 3> >();
          break;
        case itk::ImageIOBase::USHORT :
          converter=new ImageConverter<itk::VectorImage<unsigned short, 3> >();
          break;
        case itk::ImageIOBase::SHORT :
          converter=new ImageConverter<itk::VectorImage<short, 3> >();
          break;
        case itk::ImageIOBase::INT :
          converter=new ImageConverter<itk::VectorImage<int, 3> >();
          break; 
        case itk::ImageIOBase::UINT:
          converter=new ImageConverter<itk::VectorImage<unsigned int, 3> >();
           break; 
        case itk::ImageIOBase::FLOAT :
          converter=new ImageConverter<itk::VectorImage<float, 3> >();
          break; 
        case itk::ImageIOBase::DOUBLE:
          converter=new ImageConverter<itk::VectorImage<double, 3> >();
          break; 
        default:
          itk::ExceptionObject("Unsupported component type");
      }
    } else if(nd==4  ) { //|| (nd==3 && nc>1)
      if(verbose) std::cout<<"Writing 4D image..."<<std::endl;
        switch(ct)
        {
          case itk::ImageIOBase::UCHAR:
            converter=new ImageConverter<itk::VectorImage<unsigned char, 4> >();
            break;
          case itk::ImageIOBase::CHAR:
            converter=new ImageConverter<itk::VectorImage<char, 4> >();
            break;
          case itk::ImageIOBase::USHORT:
            converter=new ImageConverter<itk::VectorImage<unsigned short, 4> >();
            break;
          case itk::ImageIOBase::SHORT:
            converter=new ImageConverter<itk::VectorImage<short, 4> >();
            break;
          case itk::ImageIOBase::INT:
            converter=new ImageConverter<itk::VectorImage<int, 4> >();
            break; 
          case itk::ImageIOBase::UINT:
            converter=new ImageConverter<itk::VectorImage<unsigned int, 4> >();
            break; 
          case itk::ImageIOBase::FLOAT:
            converter=new ImageConverter<itk::VectorImage<float, 4> >();
            break; 
          case itk::ImageIOBase::DOUBLE:
            converter=new ImageConverter<itk::VectorImage<double, 4> >();
            break; 
          default:
            itk::ExceptionObject("Unsupported component type");
        }
    } else {
      throw itk::ExceptionObject("Unsupported number of dimensions");
    }

    if(converter)
    {
      if(minc_float)
          minc_type=typeid(float).name();
      else if(minc_double)
          minc_type=typeid(double).name();
      else if(minc_byte)
          minc_type=typeid(unsigned char).name();
      else if(minc_short)
          minc_type=typeid(unsigned short).name();
      else if(minc_int)
          minc_type=typeid(int).name();
      
      converter->setup(verbose,assume_dti,use_b_matrix,dwi_flip_z,inv_x,inv_y,inv_z,center,show_meta,history,minc_type,minc_to_nrrd,nrrd_to_minc);
      converter->load_and_save_image(io,output.c_str(),oct);
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
