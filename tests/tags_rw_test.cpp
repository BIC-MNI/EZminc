/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2006 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include "minc_helpers.h"
#include <iostream>
#include <getopt.h>


using namespace  minc;
using namespace  std;

void show_usage(const char *name)
{
  std::cerr 
      << "This tests tags reading and writing"<<std::endl
      << std::endl
      << "Usage: "<<name<<" <input tag> <output tag>" << std::endl
      << "--vol2 try to read second volume"<< std::endl
      << "--clobber overwrite output file" << std::endl;
}

int main (int argc, char **argv)
{
  int verbose=0;
  int clobber=0;
  int vol=1;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose,      1},
    {"clobber", no_argument, &clobber,      1},
    {"vol2", no_argument, &vol,     2 },
    {"quiet",   no_argument, &verbose,      0},
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
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) !=2)
  {
    show_usage(argv[0]);
    return 1;
  }

  std::string output=argv[argc-1]; //last argument is output file... maybe we should make it a parameter instead?
  
  if (!clobber && !access (output.c_str(), F_OK))
  {
    std::cerr << output.c_str () << " Exists!" << std::endl;
    return 1;
  }  
    
  try
  {
		tag_points tags;
    std::vector<double> labels;
    
    tag_points test_tags;
    std::vector<double> test_labels;
        
    read_tags(tags,labels,argv[optind],vol);
    
    for(int j=0;j<tags.size();j++)
    {
      test_tags.push_back(tags[j]);
      test_labels.push_back(labels[j]);
      std::cout<<tags[j][0]<<" "<<tags[j][1]<<" "<<tags[j][2]<<" "<<labels[j]<<std::endl;
    }
    
		//read_tags(tags,argv[1]);
		//read_tags(tags2,argv[2]);
		//		write_2tags(tags,tags2,argv[3]);
    
    write_tags(test_tags,test_labels,output.c_str());
	} catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
		cerr << err.msg()<<endl;
    return 1;
  }
	return 0;
}
