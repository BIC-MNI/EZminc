/* ----------------------------- MNI Header -----------------------------------
@NAME       :  volume_gtc_similarity
@DESCRIPTION:  an example of calculating generalized volume similarity metrics 
@COPYRIGHT  :
              Copyright 2011 Vladimir Fonov, McConnell Brain Imaging Centre, 
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
#include "minc_1_simple.h"
#include <getopt.h>
#include <algorithm>
#include <string.h>
//#include <map>
#include <set>

using namespace minc;

void show_usage (const char * prog)
{
  std::cout<<"Program calculates multiple volume similarity metrics for discrete labels "<<std::endl
           <<"or Generalized Tanimoto coefficient (GTC)" <<std::endl
           <<"based on :  William R. Crum, Oscar Camara, and Derek L. G. Hill"
           <<"\"Generalized Overlap Measures for Evaluation and Validation in Medical Image Analysis \""
           <<" IEEE TRANSACTIONS ON MEDICAL IMAGING, VOL. 25, NO. 11, NOVEMBER 2006"<<std::endl
           <<"http://dx.doi.org/10.1109/TMI.2006.880587"<<std::endl<<std::endl
           <<"Usage: "<<prog<<" <input1.mnc> <input2.mnc> "<<std::endl
           <<"\t[--background] take background voxels (0) into account" <<std::endl
           <<"\t[--gkappa] calculate generalized kappa "<<std::endl
           <<"\t[--gtc]    calculate generalized tinamoto overlap "<<std::endl
           <<"\t[--akappa] calculate average kappa"<<std::endl
           <<"\t[--csv]    output overla metrics in the format gkappa,gtc,akappa"<<std::endl
           <<"\t[--exclude l1[,l2[,l3]]]  exclude labels from considerations "<<std::endl
           <<"\t[--include l1[,l2[,l3]]]  consider only listed labels"<<std::endl;
}

typedef unsigned short voxel_type;

int main(int argc,char **argv)
{
  int verbose=0;
  int csv=0;
  int gkappa=0;
  int gtc=0;
  int akappa=0;
  int background=0;
  
  static struct option long_options[] = {
    {"verbose",   no_argument,            &verbose, 1},
    {"quiet",     no_argument,            &verbose, 0},
    {"gkappa",    no_argument,            &gkappa, 1},
    {"akappa",    no_argument,            &akappa, 1},
    {"gtc",       no_argument,            &gtc, 1},
    {"csv",       no_argument,            &csv, 1},
    {"bg",        no_argument,            &background, 1},
    {"background",no_argument,            &background, 1},
    {"exclude",   required_argument,      0,      'e'},
    {"include",   required_argument,      0,      'i'},
   
    {0, 0, 0, 0}
  };
  
  std::set<voxel_type> exclude;
  std::set<voxel_type> include;
  
  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vqe:i:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c)
    {
      case 0:
        break;
        
      case 'v':
        std::cout << "Version: 0.1" << std::endl;
        return 0;
        
      case 'e':
        {
          const char* delim=", ";
          
          for(char *tok=strtok(optarg,delim);tok;tok=strtok(NULL,delim))
            exclude.insert(static_cast<voxel_type>(atoi(tok)));
        }
        break;
        
      case 'i':
        {
          const char* delim=", ";
          
          for(char *tok=strtok(optarg,delim);tok;tok=strtok(NULL,delim))
            include.insert(static_cast<voxel_type>(atoi(tok)));
        }
        break;
        
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
    }
  }

  if((argc - optind) < 2) {
    show_usage (argv[0]);
    return 1;
  }
  
  if(! (gkappa||gtc||akappa) ) //no options selected, do all 
  {
    verbose=1;
    gkappa=1;
    gtc=1;
    akappa=1;
  }
  try
  {
    minc_1_reader rdr1;
    rdr1.open(argv[optind]);
    
    minc_1_reader rdr2;
    rdr2.open(argv[optind+1]);
    
    if(rdr1.dim_no()!=rdr2.dim_no() )
    {
      std::cerr<<"Different number of dimensions!"<<std::endl;
      return 1;
    }
    
    size_t size=1;
    
    for(int i=0;i<5;i++)
    {
      if(rdr1.ndim(i)!=rdr2.ndim(i))
        std::cerr<<"Different dimensions length! "<<std::endl;
      
      if(rdr1.ndim(i)>0) size*=rdr1.ndim(i);
    }
    
    for(int i=0;i<5;i++)
    {
      if(rdr1.nspacing(i)!=rdr2.nspacing(i) )
        std::cerr<<"Different step size! "<<std::endl;
    }
    
    rdr1.setup_read_ushort();
    rdr2.setup_read_ushort();
    
    std::vector<voxel_type> buffer1(size),buffer2(size);
    
    load_standard_volume<voxel_type>(rdr1,&buffer1[0]);
    load_standard_volume<voxel_type>(rdr2,&buffer2[0]);
    
    // Let's find all the unique labels in both volumes!
    std::set<voxel_type> label_set;
    if(include.empty())
    {
      for(size_t i=0; i<size ; i++ )
      {
        voxel_type vx=buffer1[i];
        if(!background && vx==0)
          continue;

        if( exclude.find(vx) != exclude.end() )
          continue; //skip this label

        label_set.insert(vx);
      }
    } else {
      // using only requested labels
      label_set=include;
    }
    
    if( verbose && !csv)
    {
      std::cout<<"Number of unique labels:"<<label_set.size()<<std::endl;
      if(background && include.empty())
        std::cout<<"Including background!"<<std::endl;
    }

    double intersect=  0.0;
    double overlap=    0.0;
    double volume=     0.0;
    double _akappa=    0.0;
    double _akappa_cnt=0.0;
    
    //tools for calculating average kappa (MICCAI2012 MultiAtlas segmentation style)
    double count_intersect,kappa;
    double count_ref,count_res;

    //TODO: add mask ? or different weighting
    for(std::set<voxel_type>::iterator it=label_set.begin();it!=label_set.end();++it)
    {
      voxel_type label=*it;
  
      count_intersect=0.0;
      count_ref=0.0;
      count_res=0.0;

      for(size_t i=0; i<size ; i++ )
      {
        if( buffer1[i]==label ) {volume+=1.0;count_ref+=1.0;}
        if( buffer2[i]==label ) {volume+=1.0;count_res+=1.0;}

        if( buffer1[i]==label && buffer2[i]==label ) {intersect+=1.0;count_intersect+=1.0;}
        if( buffer1[i]==label || buffer2[i]==label ) overlap+=1.0;
      }
      if( (count_ref+count_res)>0.0)
      {
        kappa=2.0*count_intersect/(count_ref+count_res);
        _akappa+=kappa;
        _akappa_cnt+=1.0;
      }
    }
    
    double _gkappa=2*intersect/volume;
    double _gtc=intersect/overlap;
    _akappa=_akappa/_akappa_cnt;
    
    std::cout.precision(10);

    if( csv )
    {
      std::cout<<_gkappa<<",";
      std::cout<<_gtc<<",";
      std::cout<<_akappa<<std::endl;
      
    } else {
      if( gkappa ){
        if(verbose) std::cout<<"Generalized Kappa ";
        std::cout<<_gkappa<<std::endl;
      }
      if( gtc )
      {
        if(verbose) std::cout<<"Generalized Tinamoto Coeffecient ";
        std::cout<<_gtc<<std::endl;
      }
      if( akappa )
      {
        if(verbose) std::cout<<"Average kappa ";
        std::cout<<_akappa<<std::endl;
        //std::cout<<"Count:"<<_akappa_cnt<<std::endl;
      }
    }
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  return 0;
}
