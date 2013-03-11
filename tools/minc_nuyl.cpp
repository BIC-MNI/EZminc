#include <iostream>
#include <getopt.h>
#include <unistd.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include "minc_histograms.h"
#include <time_stamp.h>    // for creating minc style history entry
#include "strtok.h"

using namespace minc;


typedef minc::simple_commulative_histogram<float> simple_histogram;

void show_usage(const char *name)
{
  std::cerr 
      << "This program implements intensity normalization algorithm published in "<<std::endl
      << "Nyul, L.G.; Udupa, J.K.; Xuan Zhang, "
      << "\"New variants of a method of MRI scale standardization,\""<<std::endl
      << "Medical Imaging, IEEE Transactions on , vol.19, no.2, pp.143-150, Feb. 2000 "<<std::endl
      << "http://dx.doi.org/10.1109/42.836373 "<<std::endl
      << std::endl
      << "Usage: "<<name<<" <source.mnc> <target.mnc> [output] or <source.mnc> <output.hist> --chist" << std::endl
      << "\t--source-mask <mask1.mnc>"<<std::endl
      << "\t--target-mask <mask2.mnc>"<<std::endl
      << "\t--steps <n> number of steps (default 10)"<<std::endl
      << "\t--fix_zero_padding fix mri volumes with lots of zeros in background"<<std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--ks calculate Kolmogorov–Smirnov difference, no output is produced"<< std::endl
      << "\t--chist output or target are comulative histograms" << std::endl;
}

void populate_chist(simple_histogram& hist,std::vector<double> &chist,std::vector<double> &pct,double cut_off,int steps,bool verbose=false)
{
  double step  = (1.0 - 2.0*cut_off/100.0) / steps ;

  chist.clear();
  pct.clear();
  chist.reserve(steps+1);
  pct.reserve(steps+1);
  //provide mapping for background
  chist.push_back(hist.min());
  pct.push_back(0.0);

  for(int i=0;i<=steps;i++)
  {
    double _pct =  i * step + cut_off/100.0;
    
    double lev_s=hist.find_percentile(_pct);
    
    if(lev_s-chist[chist.size()-1]<(hist.max()-hist.min())/1.0e5)
    {
      std::cerr<<"Warning: "<<_pct*100.0<<" percentile collapses in source, skipping"<<std::endl;
    } else {
      chist.push_back(lev_s);
      pct.push_back(_pct);
    }
  }
  //provide mapping for upper range
  chist.push_back(hist.max());
  pct.push_back(1.0);
}

void save_chist(const char* output,const std::vector<double> &chist)
{
  std::ofstream out(output);
  
  for(size_t i=0;i<chist.size();i++)
  {
    out<<chist[i];
    if(i!=(chist.size()-1))
      out<<",";
  }
  
  out<<std::endl;
}

void load_chist(const char* input,std::vector<double> &chist)
{
  std::ifstream in(input);
  char tmp[1024];

  in>>tmp;
  
  chist.clear();
  stringtok_d(chist,tmp,",");
}


int main(int argc,char **argv)
{
  int verbose=0;
  int normalize=0;
  int bimodalT=0;
  int debug=0;
  int hist_bins=4000;
  int steps=10;
  int clobber=0;
  int calc_ks=0;
  int fix_zero_padding=0;
  int chist=0;
  
  double cut_off=0.01;
  
  std::string source_mask_f,target_mask_f;
  
  char *_history = time_stamp(argc, argv); 
  std::string history=_history;
  free(_history);
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose,      1},
    {"chist",   no_argument, &chist,        1},
    {"clobber", no_argument, &clobber,      1},
    {"debug",   no_argument, &debug,        1},
    {"quiet",   no_argument, &verbose,      0},
    {"fix_zero_padding",no_argument,&fix_zero_padding,1},
    {"source-mask",   required_argument, 0, 's'},
    {"target-mask",   required_argument, 0, 't'},
    //{"bins",          required_argument, 0, 'b'},
    {"steps",         required_argument, 0, 'S'},
    {"ks", no_argument, &calc_ks,      1},
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
  
  if( chist  )
  {
    if((argc - optind)<2)
    {
      show_usage(argv[0]);
      return 1;
    }
  } else {
    if ((argc - optind) < 3 && !calc_ks || (argc - optind)<2 && calc_ks  ) {
      show_usage(argv[0]);
      return 1;
    }
  }
  
  std::string input_src_f=argv[optind];
  std::string input_trg_f=argv[optind+1];
  
  std::string output_f=(argc - optind)>2?argv[optind+2]:"";
  
  if(chist) 
  {
    std::cout<<"Using comulative histograms"<<std::endl;
    
    if((argc - optind)<3) 
    {
      output_f=argv[optind+1];
      input_trg_f="";
    } else {
      output_f=argv[optind+2];
    }
  }
  
  if (!output_f.empty() && !clobber && !access (output_f.c_str(), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }  
 
  try
  {
    simple_volume<float> src,trg;
    
    minc_1_reader rdr1;
    rdr1.open(input_src_f.c_str());
    load_simple_volume<float>(rdr1,src);
    
    simple_histogram src_hist_s;
    simple_histogram trg_hist_s;
    
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
        float src_min,src_max;
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
    std::vector<double> src_levels_s;
    std::vector<double> src_pct_s;
    
    populate_chist(src_hist_s,src_levels_s,src_pct_s,cut_off,steps,verbose);
    
    if(!input_trg_f.empty())
    {
      std::vector<double> trg_levels_s;
      std::vector<double> trg_pct_s;
      if(chist) //load commulative histogram from a file
      {
        load_chist(input_trg_f.c_str(),trg_levels_s);

      } else {
        minc_1_reader rdr2;
        rdr2.open(input_trg_f.c_str());
        load_simple_volume<float>(rdr2,trg);
        
        
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

        if(calc_ks)
        {
          if(verbose) std::cout<<"Kolmogorov–Smirnov distance:";
          double significance=0.0;
          double dist=src_hist_s.ks_distance(trg_hist_s,significance);
          std::cout<<dist<<std::endl;
          //if(verbose) std::cout<<"Significance:";
          //std::cout<<significance<<std::endl;
          return 0;
        }
        
        populate_chist(trg_hist_s,trg_levels_s,trg_pct_s,cut_off,steps,verbose);  
        
      }
      
      if(trg_levels_s.size()!=src_levels_s.size())
      {
        std::cerr<<"Source and target histogram size mismatch:"<<std::endl;
        std::cerr<<"Source:"<<src_levels_s.size()<<std::endl;
        std::cerr<<"Target:"<<trg_levels_s.size()<<std::endl;
        return 1;
      }
      
      if(verbose)
      {
        for(size_t i=0;i<src_levels_s.size();i++)
        {
          std::cout<<src_levels_s[i]<<" => "<<trg_levels_s[i]<<std::endl;
        }
        std::cout<<"Recalculating intensities..."<<std::flush;
      }
        
      
      //TODO: analyze commulative histograms for collapsing levels?
    
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
    } else {
      save_chist(output_f.c_str(),src_levels_s);
    }
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
}
