#include <iostream>
#include <getopt.h>
#include <unistd.h>

#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include "minc_histograms.h"
#include <time_stamp.h>    // for creating minc style history entry

using namespace minc;

void show_usage(const char *name)
{
  std::cerr 
      << "This program converts intensities into rank values [0-100]"
      << "Usage: "<<name<<" <source> <output>" << std::endl
      << "\t--mask <mask1.mnc> use mask to estimate histogram "<<std::endl
      << "\t--fix_zero_padding fix mri volumes with lots of zeros in background"<<std::endl
      << "\t--verbose be verbose" << std::endl;
}


typedef minc::simple_commulative_histogram<float> simple_histogram;

int main(int argc,char **argv)
{
  int verbose=0;
  int normalize=0;
  int bimodalT=0;
  int debug=0;
  int steps=10;
  int clobber=0;
  int fix_zero_padding=0;
  
  std::string mask_f;
  
  char *_history = time_stamp(argc, argv); 
  std::string history=_history;
  free(_history);
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose,      1},
    {"clobber", no_argument, &clobber,      1},
    {"debug",   no_argument, &debug,        1},
    {"quiet",   no_argument, &verbose,      0},
    {"fix_zero_padding",no_argument,&fix_zero_padding,1},
    {"mask",   required_argument, 0, 'm'},
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

  if ((argc - optind) < 2)
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
    simple_volume<float> src;
    float src_min,src_max;
    
    minc_1_reader rdr1;
    rdr1.open(input_f.c_str());
    load_simple_volume<float>(rdr1,src);
    
    simple_histogram src_hist_s;
    
    if(!mask_f.empty())
    {
      minc_byte_volume src_mask;
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<mask_f.c_str()<<std::endl;
      
      rdr.open(mask_f.c_str());
      load_simple_volume(rdr,src_mask);
      
      if(src_mask.size()!=src.size())
        REPORT_ERROR("Source mask size mismatch");
      
      src_hist_s.build_histogram(src,src_mask);
      
    } else {
      
      if(fix_zero_padding) //will try to remove zero padding
      {
        volume_min_max(src,src_min,src_max);
        
        minc_byte_volume src_mask(src.dims());
        
        for(int i=0;i<src.c_buf_size();i++)
        {
          src_mask.c_buf()[i]=src.c_buf()[i]>src_min?1:0;
        }
        src_hist_s.build_histogram(src,src_mask);
         
      } else {
        src_hist_s.build_histogram(src);
      }
    }
    
    src_min=src_hist_s.min();
    src_max=src_hist_s.max();
    //src_hist.save("src.lst");
  
    std::vector<double> src_levels_s;
    
    if(verbose) 
      std::cout<<"Recalculating intensities..."<<std::flush;
    
    for(int i=0;i<src.c_buf_size();i++)
    {
      //use LUT to map the intensities
      int bin;
      double input=src.c_buf()[i];
      
      double output=src_hist_s.rank(input);
      
      src.c_buf()[i]=output;
    }
    
    if(verbose) 
      std::cout<<"Done!"<<std::endl;
    
    minc_1_writer wrt;
    
    wrt.open(output_f.c_str(),rdr1);
    wrt.append_history(history.c_str());
    
    save_simple_volume<float>(wrt,src);
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
}