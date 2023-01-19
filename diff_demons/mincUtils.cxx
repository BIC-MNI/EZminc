#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H

#include "mincUtils.h"
#include <time.h>


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
