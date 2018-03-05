#include "minc_1_rw.h"
#include <iostream>
#include <fstream>

#include "pca_utils.h"
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <getopt.h>
#include <math.h>
#include <gsl_glue.h>
#include <gsl/gsl_eigen.h>

using namespace minc;

void save_results(const string_table& tbl, 
                  const volumes &mean, 
                  gsl_double_matrix& pc, 
                  gsl_double_vector &w, 
                  const std::string &output,
                  bool normalize)
{

  string_table out_tbl(pc.size(0)+1);
  
  for(int j=0;j<(pc.size(0)+1);j++)
  {
    out_tbl[j].resize(tbl[0].size()+1);
  }
  
  for(int k=0;k<mean.size();k++)
  {
    minc_1_reader rdr;
    rdr.open(tbl[0][k].c_str(),false,true);// I need only the metadata
    
    minc_1_writer wrt;  
    char file[1024];
    sprintf(file,"%s_%d_%03d.mnc",output.c_str(),k,0); // 0 - mean
    out_tbl[0][k+1]=minc::_basename(file);
    
    out_tbl[0][0]="0.0"; //mean

    wrt.open(file,rdr.info(),2,NC_FLOAT);
    save_simple_volume<float>(wrt,mean[k]);
    minc::simple_volume<float> tmp,tmp_read;

    for(int j=0;j<pc.size(0);j++)
    {
      char file2[1024];
      sprintf(file2,"%s_%d_%03d.mnc",output.c_str(),k,j+1); // 0 - mean
      tmp.resize(mean[k].size());
      out_tbl[j+1][k+1]=minc::_basename(file2);
      char tmp2[256];
      sprintf(tmp2,"%g",w[j]);
      out_tbl[j+1][0]=tmp2;
      tmp=0.0;
      for(int i=0;i<pc.size(1);i++)
      {
        minc_1_reader rdr2;
        rdr2.open(tbl[i][k].c_str());
        load_simple_volume<float>(rdr2,tmp_read);
        tmp_read-=mean[k];
        tmp.weighted_add(tmp_read,pc.get(i,j));
      }

      if(normalize)
        tmp/=sqrt(fabs(w[j]));

      minc_1_writer wrt2; 
      wrt2.open(file2,rdr.info(),2,NC_FLOAT);
      save_simple_volume<float>(wrt2, tmp);
      std::cout<<j<<"\t"<<std::flush;
    }
    std::cout<<std::endl;
  }
  
  std::ofstream out(output.c_str());
  for(int k=0;k<out_tbl.size();k++)
  {
    for(int j=0;j<out_tbl[0].size();j++)
    {
      out<<out_tbl[k][j];
      if(j!=(out_tbl[0].size()-1))
        out<<",";
    }
    out<<std::endl;
  }
}


void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" <training list> <output_prefix>  " << std::endl
      << "\t--mask <mask.mnc> use mask for ROI" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--normalize normalize PCA to 1" << std::endl 
      << "\t--nobias do not remove bias (mean)" << std::endl 
      << "\t--threshold <f>, default 0.98 - proportion of variance explained" << std::endl;
}

int main(int argc,char **argv)
{
  int clobber=0;
  int verbose=0;
  int normalize=0;
  int nobias=0;
  int n_select=0;
  double threshold=0.98;
  
  std::string mask_f;
  static struct option long_options[] =
  {
    {"verbose",   no_argument, &verbose, 1},
    {"quiet",     no_argument, &verbose, 0},
    {"clobber",   no_argument, &clobber, 1},
    {"normalize", no_argument, &normalize, 1},
    {"nobias",    no_argument, &nobias, 1},
    {"mask"     , required_argument, 0, 'm'},
    {"threshold", required_argument, 0, 't'},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "sm:t:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        break;
      case '?':
        /* getopt_long already printed an error message. */
      case 'm':
        mask_f=optarg;
        break;
      case 't':
        threshold=atof(optarg);
        break;
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
  
  std::string train_f=argv[optind];
  std::string output=argv[optind+1];
  
  try
  {
    string_table tbl;
    read_table(train_f.c_str(),tbl);
    if(!table_verify(tbl))
      return 1;
  
    //1. calculate averages
    volumes means;
    if(nobias)
    {
        means.resize(tbl[0].size());
        for(int i=0;i<tbl[0].size();i++) //go through modalities
        {
            minc_1_reader rdr1;
            rdr1.open(tbl[0][i].c_str());
            //HACK: just for volume resizing...
            load_simple_volume<float>(rdr1,means[i]);
            means[i]=0.0;
        }
    } else {
        std::cout<<"Averaging..."<<std::endl;
        average_tables(tbl, means);
        std::cout<<"Done"<<std::endl;
    }
    minc_byte_volume mask;
    if(!mask_f.empty())
    {
      minc_1_reader rdr;
      rdr.open(mask_f.c_str());
      load_simple_volume(rdr,mask);
    }
    
    //2. calculate correlation matrix
    gsl_double_matrix cov,pc;
    gsl_double_vector w;
    std::cout<<"Calculating covariance matrix..."<<std::endl;
    
    if(mask_f.empty())
      calculate_covariance(tbl,means,       cov, 0, means.size(), false);
    else
      calculate_covariance(tbl,means, mask, cov, 0, means.size(), false);
    
    std::cout<<"Done"<<std::endl;
    //3. PCA
    std::cout<<"Calculating PCA..., threshold="<<threshold<<std::endl;
    n_select=calculate_pca(cov, pc, w, threshold);
    std::cout<<"Done"<<std::endl;
    
    //4. store results?
    std::cout<<"Writing "<<n_select<<" eigen vectors"<<std::endl;
    if(normalize)
    {
      double ss=1;
      if(mask_f.empty())
      {
        ss=means[0].size().vol();
      } else {
        ss=mask.size().vol(); 
      }
      //ss*=ss;
      w*=1.0/ss;
    }
    save_results(tbl, means, pc, w, output, normalize);
    std::cout<<"Done"<<std::endl;
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  
  return 0;
}
