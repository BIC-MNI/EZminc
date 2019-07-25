#include <iostream>
#include <getopt.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <cstdlib>
#include "noise_estimate.h"

using namespace minc;

void show_usage(const char *name)
{
  std::cerr 
      << "This program implements noise estimation algorithm published in "<<std::endl
      << "Pierrick Coupe, Jose V. Manjon, Elias Gedamu, Douglas L. Arnold,"<<std::endl
      << "Montserrat Robles, D. Louis Collins: An Object-Based Method for Rician"<<std::endl
      << "Noise Estimation in MR Images. MICCAI (1) 2009: 601-608."<<std::endl<<std::endl
      << "Usage: "<<name<<" <input> " << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--noise  output noise level (default)" << std::endl
      << "\t--snr output SNR " << std::endl
      << "\t--bins <n> histograms bins used to detect object , default 2000" << std::endl
      << "\t--mask <mask.mnc> object mask - overrides object detection based on k-means" << std::endl
      << "\t--gauss - use gaussian noise assumption, default - disabled" << std::endl;

}


int main(int argc,char **argv)
{
  int verbose=0;
  int normalize=0;
  int bimodalT=0;
  int debug=0;
  int k_means=0;
  int maxiter=10;
  int hist_bins=2000;
  int output_snr=0; 
  int output_noise=0;
  int gaussian_noise=0;
  
  std::string input_mask_f;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose,      1},
    {"debug",   no_argument, &debug,        1},
    {"quiet",   no_argument, &verbose,      0},
    {"snr",     no_argument, &output_snr,   1},
    {"noise",   no_argument, &output_noise, 1},
    {"gauss",   no_argument, &gaussian_noise, 1},
    {"bins",   required_argument, 0, 'b'},
    {"mask",   required_argument, 0, 'm'},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "m:b:", long_options, &option_index);

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
      case 'm':
        input_mask_f=optarg;
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
  
  if(!output_snr && !output_noise) output_noise=1;
  
  try
  {
    simple_volume<float> input;
    minc_byte_volume  input_mask;
    
    minc_1_reader rdr;
    rdr.open(input_f.c_str());
    load_simple_volume<float>(rdr,input);
    
    if(!input_mask_f.empty())
    {
      minc_1_reader rdr2;
      rdr2.open(input_mask_f.c_str());
      
      for(int i=0;i<3;i++)
      {
        if(rdr.ndim(i)!=rdr2.ndim(i))
        {
          std::cerr<<"Mask has Different dimensions length! "<<std::endl;
          return 1;
        }
      }
      
      for(int i=0;i<5;i++)
      {
        if(rdr.nspacing(i)!=rdr2.nspacing(i) )
        {
          std::cerr<<"Mask has Different step size! "<<std::endl;
          return 1;
        }
      }
      load_simple_volume<unsigned char>(rdr2,input_mask);
      
      
    }
    
    double mean_signal=0.0;
    
    double nsig_corr=noise_estimate(input,mean_signal,gaussian_noise,verbose,hist_bins,input_mask);
          
    if(output_noise) {
      if(verbose) std::cout<<"Noise=";
      std::cout<<nsig_corr<<std::endl;
    }
    
    if(output_snr)
    {
      if(verbose) std::cout<<"SNR=";
      std::cout<<mean_signal/nsig_corr<<std::endl;
    }
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
}

// kate: indent-mode cstyle; indent-width 2; replace-tabs on; tab-width 2
