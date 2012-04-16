/* ----------------------------- MNI Header -----------------------------------
@NAME       : sharpness_estimate
@DESCRIPTION: estimate volume sharpness, based on median value of gradients
@COPYRIGHT  :
              Copyright 2012 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */


#include <iostream>
#include <getopt.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <stdlib.h>
#include "sharpness_estimate.h"

using namespace minc;

void show_usage(const char *name)
{
  std::cerr 
      << "This program implements simple sharpness estimate "<<std::endl
      << " based on median gradient magnitude"<<std::endl
      << "Usage: "<<name<<" <input> " << std::endl
      << "\t--mask <mask> use mask " << std::endl;

}


int main(int argc,char **argv)
{
  int verbose=0;
  std::string mask_f;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose,      1},
    {"quiet",   no_argument, &verbose,      0},
    {"mask",    required_argument, 0,      'm'},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "m:", long_options, &option_index);

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
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) < 1)
  {
    show_usage(argv[0]);
    return 1;
  }
  std::string input_f=argv[optind];
  
  try
  {
    simple_volume<float> input;
    simple_volume<unsigned char> mask;

    minc_1_reader rdr;
    rdr.open(input_f.c_str());
    load_simple_volume<float>(rdr,input);
    
    if(!mask_f.empty())
    {
      minc_1_reader rdr2;
      rdr2.open(mask_f.c_str());
      
      load_simple_volume<unsigned char>(rdr2,mask);
      if(input.size()!=mask.size())
      {
        std::cerr<<"Mask size mismatch!"<<std::endl;
        return 1;
      }
    } else {
      mask.resize(input.size());
      mask=1; //use all voxels
    }
      
    double sharpness=sharpness_estimate(input,mask);
    std::cout<<sharpness<<std::endl;
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
}
