#include <iostream>
#include <getopt.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include "minc_histograms.h"
#include <time_stamp.h>    // for creating minc style history entry

using namespace minc;

void show_usage(const char *name)
{
  std::cerr 
      << "This program implements intensity normalization algorithm published in "<<std::endl
      << "Nyul, L.G.; Udupa, J.K.; Xuan Zhang, "
      << "\"New variants of a method of MRI scale standardization,\""<<std::endl
      << "Medical Imaging, IEEE Transactions on , vol.19, no.2, pp.143-150, Feb. 2000 "<<std::endl
      << "http://dx.doi.org/10.1109/42.836373 "<<std::endl
      << std::endl
      << "Usage: "<<name<<" <source> <target> <output>" << std::endl
      << "\t--source-mask <mask1.mnc>"<<std::endl
      << "\t--target-mask <mask2.mnc>"<<std::endl
      << "\t--steps <n> number of steps (default 10)"<<std::endl
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
  int hist_bins=4000;
  int steps=10;
  int clobber=0;
  int fix_zero_padding=0;
  
  double cut_off=0.01;
  
  std::string source_mask_f,target_mask_f;
  
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
    {"source-mask",   required_argument, 0, 's'},
    {"target-mask",   required_argument, 0, 't'},
    //{"bins",          required_argument, 0, 'b'},
    {"steps",         required_argument, 0, 'S'},
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
      case 'b':
        hist_bins=atoi(optarg);
        break;
      case 's':
        source_mask_f=optarg;
        break;
      case 't':
        target_mask_f=optarg;
        break;
      case 'S':
        steps=atoi(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) < 3)
  {
    show_usage(argv[0]);
    return 1;
  }
  std::string input_src_f=argv[optind];
  std::string input_trg_f=argv[optind+1];
  std::string output_f=argv[optind+2];
  
  if (!clobber && !access (output_f.c_str(), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }  
 
  try
  {
    simple_volume<float> src,trg;
    float src_min,src_max;
    float trg_min,trg_max;
    
    minc_1_reader rdr1;
    rdr1.open(input_src_f.c_str());
    load_simple_volume<float>(rdr1,src);
    
    simple_histogram src_hist_s;
    
    if(!source_mask_f.empty())
    {
      minc_byte_volume src_mask;
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<source_mask_f.c_str()<<std::endl;
      
      rdr.open(source_mask_f.c_str());
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
  
    minc_1_reader rdr2;
    rdr2.open(input_trg_f.c_str());
    load_simple_volume<float>(rdr2,trg);
    simple_histogram trg_hist_s;
    
    if(!target_mask_f.empty())
    {
      minc_byte_volume   trg_mask;
      
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<target_mask_f.c_str()<<std::endl;
      
      rdr.open(target_mask_f.c_str());
      load_simple_volume(rdr,trg_mask);
      
      if(trg_mask.size()!=trg.size())
        REPORT_ERROR("Target mask size mismatch");
      
      trg_hist_s.build_histogram(trg,trg_mask);
    } else {
      //build histograms 
      trg_hist_s.build_histogram(trg);
    }
    trg_min=trg_hist_s.min();
    trg_max=trg_hist_s.max();
    
    double step  = (1.0 - 2.0*cut_off/100.0) / steps ;
    
    std::vector<double> src_levels_s;
    std::vector<double> trg_levels_s;
    
    //provide mapping for background
    src_levels_s.push_back(src_min);
    
    trg_levels_s.push_back(trg_min);
    
    if(verbose)
      std::cout<<"[ min ] "<<src_levels_s[0]<<" => "<<trg_levels_s[0]<<std::endl;
    
    for(int i=0;i<=steps;i++)
    {
      double pct =  i * step + cut_off/100;
      
     
      double src_lev_s=src_hist_s.find_percentile(pct);
      double trg_lev_s=trg_hist_s.find_percentile(pct);
      
      if(trg_lev_s-trg_levels_s[trg_levels_s.size()-1]<(trg_max-trg_min)/100000.0)
      {
        std::cerr<<"Warning: "<<pct*100<<" percentile collapses in target, skipping"<<std::endl;
      } else {
        
        src_levels_s.push_back(src_lev_s);
        trg_levels_s.push_back(trg_lev_s);
        
        if(verbose)
          std::cout<<"[ "<<pct*100.0<<" ] "<<src_levels_s[src_levels_s.size()-1]<<" => "<<trg_levels_s[trg_levels_s.size()-1]<<std::endl;
      }
    }
    //provide mapping for upper range
    
    src_levels_s.push_back(src_max);
    trg_levels_s.push_back(trg_max);
    
    if(verbose)
      std::cout<<"[ max ] "<<src_levels_s[src_levels_s.size()-1]<<" => "<<trg_levels_s[trg_levels_s.size()-1]<<std::endl;
    
    if(verbose) 
      std::cout<<"Recalculating intensities..."<<std::flush;
    
    for(int i=0;i<src.c_buf_size();i++)
    {
      //use LUT to map the intensities
      int bin;
      double input=src.c_buf()[i];
      
      double output=input;
      for(bin=0;bin<src_levels_s.size();bin++)
        if(input <= src_levels_s[bin]) break;
      
      if(bin==0) // first bin ?
        output=trg_levels_s[0];
      else if(bin>=(src_levels_s.size()-1))  
        output=trg_levels_s[trg_levels_s.size()-1];
      else 
        output=(input-src_levels_s[bin-1])/(src_levels_s[bin]-src_levels_s[bin-1])*(trg_levels_s[bin]-trg_levels_s[bin-1])+trg_levels_s[bin-1];
      
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