/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2009 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include "minc_1_rw.h"
#include <iostream>
#include <fstream>
#include <algorithm>

#include "pca_utils.h"
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <getopt.h>
#include <cmath>
#include <gsl_glue.h>
#include <gsl_gauss.h>

using namespace minc;

void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" <file> <lm_file 1> <lm_file 2> .. <lm_file n>  " << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask.mnc>"<<std::endl
      << "\t--fixed <n> number of fixed components, with weight 1"<<std::endl;
}

int main(int argc,char **argv)
{
  int clobber=0;
  int verbose=0;
  int normalize=0;
  std::string mask_f;
  int fixed_comp=0;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"quiet",   no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"mask",    required_argument, 0, 'm'},
    {"fixed",    required_argument, 0, 'f'},
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
      case 'f':
        fixed_comp=atoi(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) < 2)
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string input_f=argv[optind];
  std::vector<std::string> output;
 
  try
  {
    //check if output exists
    
    int ninput= argc - optind - 1;
    
    std::vector<minc_1_reader> files;
    std::vector<minc_input_iterator<float> > in;
    
    minc_1_reader mask;
    minc_input_iterator<unsigned char> mask_it;
      
    if(!mask_f.empty())
    {
      if(verbose)
        std::cout<<"Opening mask:"<<mask_f.c_str()<<std::endl;
      
      mask.open(mask_f.c_str());
      mask.setup_read_byte();
      mask_it.attach(mask);
      mask_it.begin();
    }
    
    files.resize(ninput+1);
    in.resize(ninput+1);
      
    for(int i=0;i<=ninput;i++) //open all input files
    {
      if(verbose)
        std::cout<<"Opening volume:"<<argv[optind+i]<<std::endl;
        
      files[i].open(argv[optind+i]);
      
      if(i==0 && !mask_f.empty())
      {
        if(files[i].dim_no()!=mask.dim_no())
        {
          std::cerr<<"Input file "<< argv[optind+i] <<" should have same number of dimensions as mask!"<<std::endl;
          return 1;
        }
        bool good=true;
        for(int j=0;j<files[i].dim_no();j++)
          if(mask.dim(j).length!=files[i].dim(j).length)
            good=false;
        if(!good)
        {
          std::cerr<<"Input file "<< argv[optind+i] <<" should have same dimensions as mask!"<<std::endl;
          return 1;
        }
      }
      //check to make sure that all files are proper
      if(i>0)
      {
        if(files[0].dim_no()!=files[i].dim_no())
        {
          std::cerr<<"Input file "<< argv[optind+i] <<" should have same number of dimensions as first file!"<<std::endl;
          return 1;
        }
        bool good=true;
        for(int j=0;j<files[0].dim_no();j++)
          if(files[i].dim(j).length!=files[0].dim(j).length)
            good=false;
        if(!good)
        {
          std::cerr<<"Input file "<< argv[optind+i] <<" should have same dimensions as first file!"<<std::endl;
          return 1;
        }
      }
      files[i].setup_read_float();
      in[i].attach(files[i]);
      in[i].begin();
    }
    
    int len=files[0].ndim(1)*files[0].ndim(2)*files[0].ndim(3)*(files[0].ndim(0)?files[0].ndim(0):1);
    
    int progress=0;
    int progress_steps=1;
    std::vector<int> report(100/progress_steps,0);
  

    LSQ_Gauss<double> lsq(ninput-fixed_comp); // linear model
    std::vector<double> basis(ninput-fixed_comp);
    if(verbose)
    {
      std::cout<<"Using :"<<fixed_comp<<" fixed components"<<std::endl;
    }
    
    while(!in[0].last())
    {
      bool good=true;
      //voxels.clear();
      //selected_ages.clear();
      
      if(verbose)
      {
        progress++;
        int pct=(progress*100)/len;
        if(!report[pct/progress_steps] )
        {
          std::cout<<pct<<"% "<<std::flush;
          report[pct/progress_steps]=1;
        }
      }
      
      if(!mask_f.empty())
      {
        if(mask_it.last()) break;
        good= (mask_it.value()!=0);
        mask_it.next();
      }
      
      double RMS=0.0;
      double iter=0.0;
      int cnt=0;
      
      double min=10;
      double max=0;
      
      double value=in[0].value();
      in[0].next();
      double bias=0.0;
      
      for(int i=0;i<fixed_comp;i++)
      {
        bias+=in[i+1].value();
        in[i+1].next();
      }
      
      for(int i=fixed_comp;i<ninput;i++)
      { 
        basis[i-fixed_comp]=in[i+1].value();
        in[i+1].next();
      }
      
      if(good)
      {
        lsq.accumulate(basis,value-bias);
      } else {
      }
    }
    
    if(verbose) {
      std::cout<<"Solution:" << std::endl;
    }
    
    std::vector<double> sol(ninput);
    lsq.solve(sol);
    for(int i=0;i<fixed_comp;i++)
      std::cout<<"1,";
    
    for(int i=fixed_comp;i<ninput;i++)
    {
      std::cout<<sol[i-fixed_comp];
      if(i<(ninput-1)) std::cout<<",";
    }
    
    std::cout<<std::endl;
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  
  return 0;
}
