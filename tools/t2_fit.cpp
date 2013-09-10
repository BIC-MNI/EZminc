#include <iostream>
#include <minc_1_simple.h>

#include <math.h>
#include "gsl_glue.h"
#include "gsl_gauss.h"

// for get_opt_long
#include <getopt.h>


using namespace  std;
using namespace  minc;


void show_usage(const char *name)
{
  std::cerr 
	  << "Usage: "<<name<<" <input1> .... <inputn>  <output> " << std::endl
    << "\tn should be more than 1"<< std::endl
    << "Optional parameters:" << std::endl
    << "\t--verbose be verbose" << std::endl
    << "\t--clobber clobber the output files" << std::endl
    << "\t--version print version" << std::endl 
    << "\t--threshold <f> threshold for background"<< std::endl
    << "\t--t2_threshold <f> threshold for maximum T2"<< std::endl
    << "\t--mask <minc_file> mask file"<< std::endl;
}

int main (int argc, char **argv)
{
  int clobber=0;
  int f_stat=0;
  double f_coeff=1.0;
  double threshold=0.0;
  double t2_threshold=1.0;
  int verbose=0;
  std::string mask_f;
  // read the arguments
	static struct option long_options[] =
	  {
		  {"verbose", no_argument, &verbose, 1},
		  {"quiet", no_argument, &verbose, 0},
		  {"clobber", no_argument, &clobber, 1},
		  {"threshold", required_argument, 0, 't'},
		  {"t2_threshold", required_argument, 0, 'h'},
      {"mask", required_argument, 0, 'm'},
		  {"version", no_argument, 0, 'v'},
		  {0, 0, 0, 0}
	  };

	int c;
	for (;;)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "v", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
		case 'v':
			cout << "Version: 1.0" << endl;
			return 0;
    case 't':
      threshold=atof(optarg);
      break;
    case 'h':
      t2_threshold=atof(optarg);
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

	if ((argc - optind) < 3)
	{
		show_usage(argv[0]);
		return 1;
	}
  std::string output=argv[argc-1]; //last argument is output file... maybe we should make it a parameter instead?
  argc-=optind+1;
  double max_te=0;
  
  if (!clobber && !access (output.c_str(), F_OK))
  {
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
  try
  {
    std::vector<minc_1_reader> files;
    std::vector<minc_input_iterator<float> > in;
    std::vector<double> te;
    files.resize(argc);
    in.resize(argc);
    minc_1_reader mask;
    minc_input_iterator<unsigned char> mask_it;
    
    if(!mask_f.empty())
    {
      if(verbose)
        std::cout<<"Opening "<<mask_f.c_str()<<std::endl;
      mask.open(mask_f.c_str(),false);
      mask.setup_read_byte();
      mask_it.attach(mask);
      mask_it.begin();
    }
    
    
    if(verbose) 
      std::cout<<"Opening "<<argc<<" files ..."<<std::endl;
    int i;
    te.resize(argc);
    for(i=0;i<argc;i++)
    {
      if(verbose)
        std::cout<<"Opening "<<argv[i+optind]<<std::endl;
      files[i].open(argv[i+optind],false);
      if(i==0 && !mask_f.empty())
      {
        if(files[0].dim_no()!=mask.dim_no())
        {
          std::cerr<<"Input file "<< argv[i+optind] <<" should have same number of dimensions as mask!"<<std::endl;
          return 1;
        }
        bool good=true;
        for(int j=0;j<files[0].dim_no();j++)
          if(mask.dim(j).length!=files[0].dim(j).length)
            good=false;
        if(!good)
        {
          std::cerr<<"Input file "<< argv[i+optind] <<" should have same dimensions as mask!"<<std::endl;
          return 1;
        }
      }
            
      //check to make sure that all files are proper
      if(i>0)
      {
        if(files[0].dim_no()!=files[i].dim_no())
        {
          std::cerr<<"Input file "<< argv[i+optind] <<" should have same number of dimensions as first file!"<<std::endl;
          return 1;
        }
        bool good=true;
        for(int j=0;j<files[0].dim_no();j++)
          if(files[i].dim(j).length!=files[0].dim(j).length)
            good=false;
        if(!good)
        {
          std::cerr<<"Input file "<< argv[i+optind] <<" should have same dimensions as first file!"<<std::endl;
          return 1;
        }
      }
      te[i]=files[i].att_value_double("acquisition","echo_time")[0];
      files[i].setup_read_float();
      in[i].attach(files[i]);
      in[i].begin();
      if(verbose)
        std::cout<<"Echo time:"<<te[i]<<std::endl;
      if(te[i]>max_te) max_te=te[i];
    }
    
    minc_info output_info=files[0].info();
    if(verbose)
      std::cout<<"Writing to:"<<output.c_str()<<std::endl;
    minc_1_writer wrt;
    wrt.open(output.c_str(),output_info,2,NC_FLOAT);
    wrt.setup_write_float();
    minc_output_iterator<float> out(wrt);
    
    if(verbose) 
      std::cout<<"Calculating T2 fitting..."<<std::endl<<std::flush;

    // total number of elements, I don't care about coordinates here
    int len=files[0].ndim(1)*files[0].ndim(2)*files[0].ndim(3);
    
    int count_good=0;
    std::vector<double> vec(argc);
    MNK_Gauss_Polinomial fit(2); // linear approximation y=a+b*x
    std::vector<double> sol(2);
    double max_t2=-100;
    for(out.begin();!out.last();out.next())
    {
      bool good=true;
      if(!mask_f.empty())
      {
        if(mask_it.last()) break;
        good= (mask_it.value()!=0);
        mask_it.next();
      }
      fit.clear();
      double iT2=0.0;
      double T2=0.0;
      for(int i=0;i<argc;i++)
      { 
        double v=in[i].value();
        
        if(v<=threshold || !good)
          good=false;
        else
          fit.accumulate(-te[i],log(v)); 
        in[i].next();
      }
      
      if(good) 
      {
        fit.solve(sol);
        iT2=sol[1];
        
        if(iT2> (1.0/t2_threshold) )
        {
          T2=1.0/iT2;
          count_good++;
          if(T2>max_t2) max_t2=T2;
        } else {
          T2=max_t2>0?max_t2:max_te;
        }
        
      }
      out.value(T2);
    }
    if(verbose) {
      std::cout<<"done!" << std::endl;
      std::cout<<count_good*100.0/len<<"% with meaningfull results" << std::endl;
    }
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg() << std::endl;
    return 1;
  }
  return 0;
}
