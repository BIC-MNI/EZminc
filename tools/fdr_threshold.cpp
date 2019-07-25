#include <iostream>
#include <minc_1_rw.h>
#include <minc_1_simple.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple_rw.h>
#include <getopt.h>
#include <cstdlib>

#include <cmath>
//#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <algorithm>

using namespace std;
using namespace minc;
void print_ref(void);

void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<"  <input> <df> <significance %%> [output]" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask>" <<std::endl
      << "\t--twosided perform two sided test, otherwise only check positive t"<<std::endl
      << "\t--negative invert t for single-sided test"<<std::endl
      << "\t--pval print corrected P value instead of t"<<std::endl
      << "\t--fstats use F-statistics instead of t"<<std::endl
      << "\t--df2 <n> 2nd degree of freedom for F statistics, by default assumes df-1"<<std::endl
      << "\t--independent assume independent samples  (less conservative test)"<<std::endl
      << "\t--out <output-pval.mnc> output p-values"<<std::endl
      << std::endl<<"Warning: by default more conservative (non dependent samples) assumption is used!"<<std::endl<<std::endl;
  print_ref();
}

void print_ref(void)
{
  std::cout<<"Christopher R. Genovese, Nicole A. Lazar, Thomas Nichols, Thresholding of Statistical Maps in Functional Neuroimaging Using the False Discovery Rate, NeuroImage, Volume 15, Issue 4, April 2002, Pages 870-878, ISSN 1053-8119, DOI: 10.1006/nimg.2001.1037."<<std::endl;
}

struct pt
{
  double p;
  double t;
  pt(double _p,double _t):p(_p),t(_t)
  {
  }
};

class ascending_sort {
  public:
    bool operator()(const pt &i,const pt &j) const { return (i.p<j.p);}
};

int main ( int argc, char **argv )
{
  std::string mask_f;
  std::string out_f;
  
  int verbose=0;
  int clobber=0;
  int independent=0;
  int twosided=0;
  int sign=1;
  int pval=0;
  int fstats=0;
  int df=0,df1=0,df2=0;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument,   &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"independent", no_argument, &independent, 1},
    {"fstats", no_argument, &fstats, 1},
    {"mask", required_argument, 0, 'm'},
    {"out", required_argument, 0, 'o'},
    {"df2", required_argument, 0, 'd'},
    {"twosided", no_argument, &twosided, 1},
    {"negative", no_argument, &sign, -1},
    {"pval", no_argument, &pval, 1},
    {0, 0, 0, 0}
  };
  
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "m:", long_options, &option_index);

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
      case 'o':
        out_f=optarg;
        break;
      case 'd':
        df2=atoi(optarg);
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
  std::string input_f=argv[optind];
  std::string output_f;
  
  df=df1=atoi ( argv[optind+1] );
  if(df2<=0) df2=df-1;
  
  double perc=atof ( argv[optind+2] )/100.0;
  
  if( (argc-optind)>3)
    output_f=argv[optind+3];
  
	try
	{

    int i,j,k;

		minc_float_volume vol;
		minc_1_reader rdr1;
		rdr1.open ( input_f.c_str() );
    
		load_simple_volume<float> ( rdr1,vol );
    
    
//put voxelvalues in a vector
    vector <pt> tval;
    minc_byte_volume mask(vol.size());
    mask=1;
    
    
    if(mask_f.empty()) //use the whole volume
    {
      
      for ( i=0; i<vol.c_buf_size(); i++ )
      {
        double t;
        if(twosided)
          t=fabs(vol.c_buf()[i]);
        else
          t=sign*vol.c_buf()[i];
          
        tval.push_back( pt(1.0 - gsl_cdf_tdist_P(t,df),t) );
      }
    } else {
      minc_1_reader rdr2;
      rdr2.open ( mask_f.c_str() );
      
      if(rdr1.dim_no()!=rdr2.dim_no() )
      {
        std::cerr<<"Different number of dimensions!"<<std::endl;
        return 1;
      }
      unsigned long size=1;
    
      for(i=0;i<5;i++)
      {
        if(rdr1.ndim(i)!=rdr2.ndim(i))
          std::cerr<<"Different dimensions length! "<<std::endl;
      
        if(rdr1.ndim(i)>0) size*=rdr1.ndim(i);
      }
    
      for(i=0;i<5;i++)
      {
        if(rdr1.nspacing(i)!=rdr2.nspacing(i) )
          std::cerr<<"Different step size! "<<std::endl;
      }
    
      load_simple_volume<unsigned char> ( rdr2,mask );
      
      for ( i=0; i<vol.c_buf_size(); i++ )
      {
        if(mask.c_buf()[i]>0)
        {
          double t;
          if(twosided)
            t=fabs(vol.c_buf()[i]);
          else
            t=sign*vol.c_buf()[i];
          
          if(fstats)
            tval.push_back( pt(1.0 - gsl_cdf_fdist_P(t,df1-df2,df2),t) );
          else
            tval.push_back( pt(1.0 - gsl_cdf_tdist_P(t,df),t) );
        }
      }
    }
    
//change voxelvalues (t-statistics) into p-values
    if(verbose)
    {
      std::cout<<"Calculating FDR threshold:"<<std::endl;
      print_ref();
      
      
      std::cout<<"Number of voxels:"<<tval.size()<<std::endl;
      std::cout<<"Df="<<df<<std::endl;
      
      if(fstats)
        std::cout<<"Warning: using F-statistics instead of T-statistics!"<<std::endl;
    }
    
    std::vector<double> fcumdist(tval.size());
    double c=0.0;
    
    for ( k=0; k<tval.size(); k++ )
      c+=1.0/(1+k);
    
    if(!out_f.empty())
    {
      minc_float_volume pvol(vol.size());
      
      int j=0;
      for ( i=0; i<vol.c_buf_size(); i++ )
      {
        if(mask.c_buf()[i]>0)
          pvol.c_buf()[i]=tval[j++].p;
        else
          pvol.c_buf()[i]=0.0;
      }
      
      if(verbose)
        std::cout<<"Saving probabilitys to :"<<out_f<<std::endl;
      
      minc_1_writer wrt1;
      wrt1.open ( out_f.c_str(),rdr1.info(),2,NC_FLOAT );

      save_simple_volume<float> ( wrt1,pvol );
    }
    
    //sort p-values from low to high to get 'cumulative distribution'
    //ascending_sort asc;
    std::sort ( tval.begin(),     tval.end(),ascending_sort() );
    
    if(independent) {
      if(verbose)
        std::cout<<"Using independent samples assumption"<<std::endl;
      c=1;
    } else {
      if(verbose)
        std::cout<<"Using non-independent samples assumption (more conservative)"<<std::endl;
    }
    
    if(twosided) {
      if(verbose)
        std::cout<<"Using two sided test, adjusting p "<<std::endl;
      perc/=2.0;
    }
    
    double slope=perc/c;
    
    if(verbose)
    {
      std::cout<<"Lowest p=" <<tval[0].p<<std::endl;
      std::cout<<"Highest p="<<tval[tval.size()-1].p<<std::endl;
      
      if(twosided)
      {
        std::cout<<"Largest absolute t=" <<tval[0].t<<std::endl;
        std::cout<<"Smallest absolute t="<<tval[tval.size()-1].t<<std::endl;
      } else {
        std::cout<<"Largest "<<(fstats?"F=":"t=") <<tval[0].t*sign<<std::endl;
        std::cout<<"Smallest "<<(fstats?"F=":"t=")<<tval[tval.size()-1].t*sign<<std::endl;
      }
      
      std::cout<<"q/c="<<slope<<std::endl;
    }
      
    //find the first p value which is above the threshold
    for (k=0; k<tval.size(); k++ )
    {
      if(tval[k].p>((double)(k+1)*slope/tval.size()))
        break;
    }
    
    double threshold_p=0.0;
    double threshold_t=10.0;
    
    if(k>0)
    {
      threshold_p=tval[k-1].p;
      threshold_t=tval[k-1].t;
      
      if(verbose)
      {
        std::cout<<"Significance threshold="<<threshold_p*100<<"%"<<std::endl;
        std::cout<<"Correspondig "<<(fstats?"F=":"t=")<< threshold_t<<std::endl;
      } else {
        if(pval)
          std::cout<<threshold_p*100<<"%"<<std::endl;
        else
          std::cout<<threshold_t<<std::endl;
      }
    } else {
      if(verbose)
        std::cout<<"Nothing is significant"<<std::endl;
      else 
        std::cout<<"NaN"<<std::endl;
      //return 1;
    }
    //apply threshold
    if(!output_f.empty())
    {
      for ( i=0; i<vol.c_buf_size(); i++ )
      {
        if(mask.c_buf()[i]<1) continue;
        
        if(twosided)
        {
          if( fabs(vol.c_buf()[i])<threshold_t)
              mask.c_buf()[i]=0;
        } else {
          if( (vol.c_buf()[i]*sign)<threshold_t)
            mask.c_buf()[i]=0;
        }
      }
      minc_1_writer wrt1;
      wrt1.open ( output_f.c_str(),rdr1.info(),2,NC_BYTE );

      save_simple_volume<unsigned char> ( wrt1,mask );
    }
	}
	catch ( const minc::generic_error & err )
	{
		std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
		std::cerr << err.msg() <<std::endl;
		return 1;

	}

	return 0;
}
