#include "mincUtils.h"
#include <time.h>

std::string minc_history;

void mincify ( itk::Object* image, const std::string & history,const char * store_datatype,itk::Object* metadata  )
 {
  if(metadata)
    image->SetMetaDataDictionary(metadata->GetMetaDataDictionary());

  if(store_datatype)
    itk::EncapsulateMetaData<std::string>(image,"storage_data_type",minc_storage_type);
  if(!history.empty())
  {
    std::string old_history;
    itk::ExposeMetaData<std::string >(image,"history",old_history);
    old_history+="\n";
    old_history+=history;
    itk::EncapsulateMetaData<std::string>(image,"history",old_history);
  }
}

//! check if the filename is XFM file, and generate required additional file names (doesn't read xfm file)
bool parse_xfm_file_name(const std::string & fname,std::string &def_field,std::string &def_field_base)
{
    int pos = fname.rfind ( "." );
    
    if ( pos != std::string::npos && pos != 0 )
    {
      if ( fname.substr ( pos ) == ".xfm" ) //trying to produce minc xfm file
      {
        def_field = fname.substr ( 0, pos );
        def_field += "_grid_0.mnc";
	
        def_field_base=(def_field); 

        size_t pos2 = def_field.rfind ( "/" );

        if ( pos2 != std::string::npos )
          def_field_base = def_field.substr ( pos2+1 );
        return true;
      }
    }
    return false;
}


std::string minc_timestamp(char **argv,int argc)
{
  std::string timestamp;
  
  char cur_time[200];
  time_t t;
  struct tm *tmp;

  t = time(NULL);
  tmp = localtime(&t);
  
  strftime(cur_time, sizeof(cur_time), "%a, %d %b %Y %T %z>>>", tmp)
  /* Get the time, overwriting newline */
  timestamp=cur_time;
  
  /* Copy the program name and arguments */
  for (i=0; i<argc; i++) {
    timestamp+=argv[i];
    timestamp+=" ";
  }
  timestamp+="\n";
}