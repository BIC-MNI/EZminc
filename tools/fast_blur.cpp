#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <math.h>
#include <limits>
#include <string>
#include <vector>
#include <getopt.h>
#include <complex>
#include <stdlib.h>
#include <minc_1_simple.h> // simple minc reading & writing
#include <minc_io_simple_volume.h>
#include <time_stamp.h>    // for creating minc style history entry
#include "fftw_blur.h"
#include "minc_1_simple_rw.h"

using namespace  std;
using namespace  minc;


void show_usage (const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" <input> <output> " << endl
      << "--clobber clobber output files"<<endl
      << "--verbose be verbose"<<endl
      << "--fwhm <f>"<< endl
      << "--dx diffirentiate in X dir"<<endl
      << "--dy diffirentiate in Y dir"<<endl
      << "--dz diffirentiate in Z dir"<<endl
      << "--grad calculate gradient (output vector field)"<<endl
      << "--gmag calculate gradient magnitude"<<endl
      << "--float output in float format "<<endl;
}

int main (int argc, char **argv)
{
  int verbose=0;
  int clobber=0;
  int low=0,hi=0;
  int dx=0,dy=0,dz=0,gmag=0,grad=0,out_float=0;
  double fwhm=1.0;
  char *history = time_stamp(argc, argv); //maybe we should free it afterwards

  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"dx", no_argument, &dx, 1},
    {"dy", no_argument, &dy, 1},
    {"dz", no_argument, &dz, 1},
    {"gmag", no_argument, &gmag, 1},
    {"grad", no_argument, &grad, 1},
    {"float", no_argument, &out_float, 1},
    {"fwhm",  required_argument, 0, 'f'},
    {0, 0, 0, 0}
    };
    
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "f:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
      break;
    case 'f':
      fwhm = atof(optarg);
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
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
  
  try {
    minc::minc_1_reader in_volume;
    
    if(verbose) std::cout<<"Reading:"<<argv[optind]<<std::endl;
    in_volume.open(argv[optind]);
    
    nc_type out_file_type=out_float?NC_FLOAT:in_volume.datatype();
    int out_file_signed=out_float?1:in_volume.is_signed();
    
    simple_volume<float> in;
    load_simple_volume<float>(in_volume,in);
    
    
    if(verbose)
        std::cout<<"FWHM:"<<fwhm/in_volume.nspacing(1)
          <<"x"<<fwhm/in_volume.nspacing(2)<<"x"<<fwhm/in_volume.nspacing(3)<<std::endl;
    
    if(grad)
    {
      minc::minc_grid_volume out(in.size());
      calc_gradient(in,out,fwhm/in_volume.nspacing(1),
                    fwhm/in_volume.nspacing(2),fwhm/in_volume.nspacing(3));
      
      minc::minc_1_writer out_volume;
      minc::minc_info out_info=in_volume.info();
      
      out_info.push_back(minc::dim_info(3,0,0,minc::dim_info::DIM_VEC));
      
      out_volume.open(argv[optind+1],out_info,3,out_file_type,out_file_signed);
      
      out_volume.copy_headers(in_volume);
      out_volume.append_history(history);
      save_simple_volume<fixed_vec<3,float> >(out_volume,out);
        
    } else if(gmag) {
        simple_volume<float> mout(in); 
        
        calc_gradient_mag(in,mout,fwhm/in_volume.nspacing(1),
                      fwhm/in_volume.nspacing(2),fwhm/in_volume.nspacing(3));
                      
        minc::minc_1_writer out_volume;
        out_volume.open(argv[optind+1],in_volume.info(),2,out_file_type,out_file_signed);
        out_volume.copy_headers(in_volume);
        out_volume.append_history(history);
        save_simple_volume<float>(out_volume,mout);
    } else {
      simple_volume<float> out(in); 
      blur_volume(in,out,dx,dy,dz,
                  fwhm/in_volume.nspacing(1),
                  fwhm/in_volume.nspacing(2),
                  fwhm/in_volume.nspacing(3));
      minc::minc_1_writer out_volume;
      out_volume.open(argv[optind+1],in_volume.info(),2,out_file_type,out_file_signed);
      out_volume.copy_headers(in_volume);
      out_volume.append_history(history);
      save_simple_volume<float>(out_volume,out);
    }
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    cerr << "errno="<<errno<<std::endl;
    free(history);
    return 1;
  }
  free(history);
  return 0;
}
