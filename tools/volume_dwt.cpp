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
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>

#include <iostream>
#include <stdlib.h>
#include <iostream>
#include <getopt.h>

#include "dwt.h"
#include "dwt_utils.h"


typedef minc::fixed_vec<3,int> idx;


using namespace minc;

void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" --forward <input>  <output_LLL> <output_LLH> <output_LHL> <output_LHH> <output_HLL> <output_HLH> <output_HHL> <output_HHH>" << std::endl
      
      << "Or: "<<name<<" --backward <input_LLL> <input_LLH> <input_LHL> <input_LHH> <input_HLL> <input_HLH> <input_HHL> <input_HHH> <output>"<<std::endl 
      
      <<"Or: "<<name<<"[--forward|--backward] --packed <input> <output>"<<std::endl;
}


int main(int argc, char **argv)
{
  int clobber=0;
  int verbose=0;
  int debug=0;
  int forward=1;
  int packed=0;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"debug",   no_argument, &debug, 1},
    {"quiet",   no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"forward", no_argument, &forward, 1},
    {"backward", no_argument, &forward, 0},
    {"packed",   no_argument, &packed, 1},
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
      /*case 'm':
        mask_f=optarg;
        break;*/
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) < (packed?2:9))
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string input_f[8],output_f[8];
  if( packed )
  {
    input_f[0]=argv[optind];
    output_f[0]=argv[optind+1];
    if (!clobber && !access (output_f[0].c_str (), F_OK))
    {
      std::cerr<<"File "<<output_f[0].c_str()<<" exists, use --clobber"<<std::endl;
      return 1;
    }    
  } else {
    if( forward )
    {
      input_f[0]=argv[optind];
      for(int i=0;i<8;i++)
      {
        output_f[i]=argv[optind+1+i];
        if (!clobber && !access (output_f[i].c_str (), F_OK))
        {
          std::cerr<<"File "<<output_f[i].c_str()<<" exists, use --clobber"<<std::endl;
          return 1;
        }
      }
    } else {
      for(int i=0;i<8;i++)
      {
        input_f[i]=argv[optind+i];
      }
      output_f[0]=argv[optind+8];
      if (!clobber && !access (output_f[0].c_str (), F_OK))
      {
        std::cerr<<"File "<<output_f[0].c_str()<<" exists, use --clobber"<<std::endl;
        return 1;
      }
    }
  }
  
  try
  {
      simple_volume<float> input_vol;
      idx input_size;
      idx padded_size;
      idx output_size;
      idx pad;
      
      minc::minc_info new_info;
      
      if(packed)
      {
        simple_volume<float> tmp_vol;
        minc_1_reader rdr;
        rdr.open(input_f[0].c_str());
        new_info=rdr.info();
        
        load_simple_volume<float>(rdr,tmp_vol);
        input_size=tmp_vol.size();
        padded_size=find_nearest_square_pow2(input_size);
        output_size=padded_size;
        
        if(forward)
        {
          input_vol.resize(padded_size);
          minc::pad_volume<float>(tmp_vol,input_vol,0); //pad volume with zeros
          
          for(int i=1;i<4;i++)
          {
            new_info[rdr.map_space(i)].start-=(padded_size[i-1]-input_size[i-1])/2;
            new_info[rdr.map_space(i)].length+=padded_size[i-1]-input_size[i-1];
          }
          
        } else {
          if(padded_size!=input_size)
          {
            std::cerr<<"Error: your volume should have sizes of the power of 2"<<std::endl;
            return 1;
          }
        }
        volume_dwt(input_vol,forward);
        
        minc_1_writer wrt;
        wrt.open(output_f[0].c_str(),new_info,2,NC_FLOAT);
        save_simple_volume<float>(wrt,input_vol);
        
      } else {
        if(forward)
        {
          simple_volume<float> tmp_vol;
          minc_1_reader rdr;
          rdr.open(input_f[0].c_str());
          new_info=rdr.info();
          load_simple_volume<float>(rdr,tmp_vol);
          
          for(int i=1;i<4;i++)
          {
            new_info[rdr.map_space(i)].start-=new_info[rdr.map_space(i)].step/2;
            new_info[rdr.map_space(i)].step*=2;
            new_info[rdr.map_space(i)].start+=new_info[rdr.map_space(i)].step/2;    
            new_info[rdr.map_space(i)].length=(new_info[rdr.map_space(i)].length+1)/2;
          }
          std::vector<simple_volume<float> > out;
          dwt_forward(tmp_vol,out);
          for(int j=0;j<8;j++)
          {
            
            minc_1_writer wrt;
            wrt.open(output_f[j].c_str(),new_info,2,NC_FLOAT);
            save_simple_volume<float>(wrt,out[j]);
          }
        } else {
          
          std::vector<simple_volume<float> > in(8);
         
          for(int j=0;j<8;j++)
          {
            minc_1_reader rdr;
            rdr.open(input_f[j].c_str());
            load_simple_volume<float>(rdr,in[j]);
            
            if(!j)
            {
              new_info=rdr.info();
          
              for(int i=1;i<4;i++)
              {
                new_info[rdr.map_space(i)].start-=new_info[rdr.map_space(i)].step/2;
                new_info[rdr.map_space(i)].step/=2;
                new_info[rdr.map_space(i)].start+=new_info[rdr.map_space(i)].step/2;    
                new_info[rdr.map_space(i)].length*=2;
              }
            }
          }
          simple_volume<float> out;
          dwt_backward(in,out);
          
          minc_1_writer wrt;
          wrt.open(output_f[0].c_str(),new_info,2,NC_FLOAT);
          save_simple_volume<float>(wrt,out);
          
        }
      }
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  return 0;
}
