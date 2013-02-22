#ifndef __mincUtils_h__
#define __mincUtils_h__

#include <itkObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>

#include <string>

//! generate minc-style timestamp
std::string minc_timestamp(int argc,char **argv);

//! change output minc file type, add history
void mincify( itk::Object* image, const std::string & history="",const char * store_datatype=typeid(float).name(),itk::Object* metadata=NULL );

//! check if the filename is XFM file, and generate required additional file names (doesn't read xfm file)
bool parse_xfm_file_name(const std::string & fname,std::string &def_field,std::string &def_field_base);

//! a helper function for image saving
template<class T> void save_image(typename T::Pointer image,const char *fname,const char * minc_datatype=NULL,itk::Object* metadata=NULL,const std::string& history="")
{
  mincify(image,history,minc_datatype,metadata);
  
  typename itk::ImageFileWriter< T >::Pointer writer = itk::ImageFileWriter<T>::New();
  writer->SetFileName(fname);
  writer->SetInput( image );
  writer->Update();
}

//! a helper function for image loading
template <class T> typename T::Pointer load_image(const char *file)
{
  typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
  reader->SetFileName(file);
  reader->Update();
  return reader->GetOutput();
}



#endif //__mincUtils_h__
