#include <minc_1_simple.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>
//#include "data_proc.h"
#include <sys/time.h>
#include <time.h>
//#include <mt19937ar.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace  std;
using namespace  minc;

void show_usage (const char *name)
{
  cerr<<"Usage:"<<name<<" <sample.mnc> <output.mnc> "<<std::endl
      <<"will generate random (uniform distribution) volume with the same parameters as sample.mnc"<<std::endl
			<<"[--clobber] clobber the output files" << std::endl
      <<"[--min <f>] minimum value (default -1)" <<std::endl
      <<"[--max <f>] maximum value (default  1)" <<std::endl
      <<"[--mask] input volume is treated as binary mask" <<std::endl
      <<"[--binary] random volume has binary distribution (50% of min or max)" <<std::endl
      <<"[--float] store volume in float (default)"<<std::endl
      <<"[--short] store volume in short"<<std::endl
      <<"[--byte] store volume in byte"<<std::endl
      <<"[--gauss <sigma> generate gaussian distribution with zero mean and sigma standard deviation, ignore min and max"<<std::endl;
}

int main (int argc, char **argv)
{
  int clobber=0;
  int mask=0,binary=0;
	double sigma=-1.0;
  const gsl_rng_type * rng_type=gsl_rng_default;
  unsigned long int rng_seed=0;
  gsl_rng * rng=NULL;
	
	int store_float=0,store_short=0,store_byte=0;
  
	try
  {
    static struct option long_options[] = { 
		  {"clobber", no_argument, &clobber, 1},
      {"min",  required_argument,   0, 'i'},
      {"max",  required_argument,   0, 'a'},
			{"gauss",  required_argument, 0, 'g'},
      {"mask", no_argument,   &mask,   1},
      {"binary", no_argument, &binary, 1},
      {"float", no_argument, &store_float, 1},
			{"short", no_argument, &store_short, 1},
			{"byte", no_argument, &store_byte, 1},
      {0, 0, 0, 0}
    };
    float min=-1.0,max=1.0;
    for (;;) {
        /* getopt_long stores the option index here. */
        int option_index = 0;
  
        int c = getopt_long (argc, argv, "a:i:", long_options, &option_index);
  
        /* Detect the end of the options. */
        if (c == -1) break;
  
        switch (c)
        {
        case 0:
          break;
        case 'a':
          max=atof(optarg);
          break;
        case 'i':
          min=atof(optarg);
          break;
				case 'g':
					sigma=atof(optarg);
					break;
        case '?':
          /* getopt_long already printed an error message. */
        default:
          show_usage (argv[0]);
          return 1;
        }
    }
  	if ((argc - optind) < 2) {
      show_usage (argv[0]);
      return 1;
    }
    //timeval timer;
    //gettimeofday(&timer, NULL);
    //init_genrand((unsigned long)timer.tv_usec);
    
		if (!clobber && !access (argv[optind+1], F_OK))
		{
			cerr << argv[optind+1] << " Exists!" << endl;
			return 1;
		}
    
    gsl_rng_env_setup();
    
    FILE *seed=fopen("/dev/urandom","rb");
    if(!seed)
    {
      std::cerr<<"Can't open /dev/urandom !\n"<<std::endl;
      return 1;
    }
    //unsigned long val[4];
    if(fread(&rng_seed,sizeof(rng_seed),1,seed)!=1)
    {
      std::cerr<<"Can't read /dev/urandom !\n"<<std::endl;
      return 1;
    }
    fclose(seed);
    //init_by_array(val,4);
    rng = gsl_rng_alloc (rng_type);
    if(!rng)
    {
      std::cerr<<"Unable to initialize GSL random number generator!"<<std::endl;
      return 1;
    }
    gsl_rng_set(rng,rng_seed);
    minc_1_reader rdr;
    rdr.open(argv[optind],true,!mask);
    unsigned int size=1;
    for(int i=0;i<5;i++)
      size*=rdr.ndim(i)>0.0?rdr.ndim(i):1;
    std::vector<float> volume(size);
    if(!mask)
    {
      for(std::vector<float>::iterator it=volume.begin();it!=volume.end();++it) {
				if(sigma>0.0)
				{ 
					*it=gsl_ran_gaussian(rng,sigma);
				} else {
					if(!binary)
						*it=gsl_rng_uniform(rng)*(max-min)+min;
					else
						*it=gsl_rng_uniform(rng)>0.5?max:min;
				}
      }
    }else{
      std::vector<unsigned char> mask(size);
      rdr.setup_read_byte();
      load_standard_volume(rdr,&mask[0]);
      std::vector<float>::iterator it;
      std::vector<unsigned char>::iterator itm;
      for(it=volume.begin(),itm=mask.begin();it!=volume.end();++it,++itm) 
      {
        if(*itm) 
        {
					if(sigma>0.0)
					{ 
						*it=gsl_ran_gaussian(rng,sigma);
					} else {
						if(!binary)
							*it=gsl_rng_uniform(rng)*(max-min)+min;
						else
							*it=gsl_rng_uniform(rng)>0.5?max:min;
					}
        }else{
          *it=0.0;
        }
      }
    }
    minc_1_writer wrt;
    wrt.open(argv[optind+1],rdr.info(),rdr.dim_no()-1,store_byte?NC_CHAR:(store_short?NC_SHORT:NC_FLOAT),store_byte?0:1);
    wrt.setup_write_float();
    save_standard_volume(wrt,&volume[0]);
    
    gsl_rng_free (rng);

		return 0;
	} catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  }
	return 0;
}
