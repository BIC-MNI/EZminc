#ifndef __mincUtils_h__
#define __mincUtils_h__

//MINC support
#if ( ITK_VERSION_MAJOR <4 ) 
#include <itkMincImageIOFactory.h>
#include <itkMincImageIO.h>
#include <itkMincHelpers.h>
#include <itkMincGeneralTransform.h>

#include <time_stamp.h>    // for creating minc style history entry
//! change output minc file type, add history
void mincify ( itk::Object* image, nc_type datatype = NC_FLOAT );
#endif
//! hold minc history
extern std::string minc_history;


//! check if the filename is XFM file, and generate required additional file names (doesn't read xfm file)
bool parse_xfm_file_name(const std::string & fname,std::string &def_field,std::string &def_field_base);


#endif //__mincUtils_h__