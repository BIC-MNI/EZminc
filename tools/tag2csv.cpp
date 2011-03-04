/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2010 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include "minc_wrappers.h"
#include <iostream>
#include <getopt.h>

using namespace  minc;
using namespace  std;

void show_usage(const char *name)
{
  std::cerr 
      << "This program converts minc .tag file into a comma-separated file"<<std::endl
      << std::endl
      << "Usage: "<<name<<" <input> <output>" << std::endl
      << "\t--clobber overwrite output file" << std::endl
      << "\t--verbose produce verbose output" << std::endl
      << "\t--[no]labels output tag labels in 4th column" << std::endl
      << "\t--[no]header print header row"<<std::endl;
}


int main (int argc, char **argv)
{
  int clobber=0;
  int labels=1;
  int header=0;
  int verbose=0;
  
  static struct option long_options[] =
  {
    {"clobber", no_argument,   &clobber,                 1},
    {"quiet",   no_argument,   &verbose,                 0},
    {"verbose",   no_argument,   &verbose,               1},
    {"labels",  no_argument,   &labels,                  1},
    {"nolabels",no_argument,   &labels,                  0},
    {"header",  no_argument,   &header,                  1},
    {"noheader",no_argument,   &header,                  0},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "r:", long_options, &option_index);

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
        show_usage(argv[0]);
        return 1;
    }
  }
  if ((argc - optind) < 2 )
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string input_f=argv[optind];
  std::string output_f=argv[optind+1]; 
  
  if (!clobber && !access (output_f.c_str(), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }
  
  
 
  try
  {
    tag_points input;
    std::vector<double> lb;
    if(labels)
        read_tags(input,lb,input_f.c_str()); 
      else
        read_tags(input,input_f.c_str()); 
    
    if(!input.size())
    {
      std::cerr<<"Warning: empty tag set!"<<std::endl;
      return 1;
    }
    std::ofstream out(output_f.c_str(),std::ios::trunc);
    if(header)
    {
      if(labels)
        out<<"X,Y,Z,LAB"<<std::endl;
      else
        out<<"X,Y,Z"<<std::endl;
    }
    for(int i=0;i<input.size();i++)
    {
      out<<input[i][0]<<","
         <<input[i][1]<<","
         <<input[i][2];
      if(labels)
        out<<","<<lb[i]<<std::endl;
      else
        out<<std::endl;
    }
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    cerr << err.msg()<<endl;
  }
}