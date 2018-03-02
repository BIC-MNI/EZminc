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
#include <gsl/gsl_linalg.h>

using namespace minc;

void save_results(const string_table& tbl, 
                  const volumes &mean, 
                  gsl_double_matrix &pc, 
                  gsl_double_vector &w, 
                  const std::string &output)
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
    
    
    char file[1024];
    sprintf(file,"%s_%d_%03d.mnc",output.c_str(),k,0); // 0 - mean
    out_tbl[0][k+1]=minc::_basename(file);
    
    out_tbl[0][0]="0.0"; //mean
    
    {
      minc_1_writer wrt;  
      wrt.open(file,rdr.info(),2,NC_FLOAT);
      save_simple_volume<float>(wrt,mean[k]);
    }
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
      minc_1_writer wrt2; 
      wrt2.open(file2,rdr.info(),2,NC_FLOAT);
      save_simple_volume<float>(wrt2,tmp);
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
      << "\t--verbose be verbose"               << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--threshold <f> threshold for the principal components, default 0.98" <<std::endl
      ;
}

int main(int argc,char **argv)
{
  int clobber=0;
  int verbose=0;
  double threshold=0.98;
  
  static struct option long_options[] =
  {
    {"verbose",   no_argument,       &verbose,  1},
    {"quiet",     no_argument,       &verbose,  0},
    {"clobber",   no_argument,       &clobber,  1},
    {"threshold", required_argument, 0,       'r'},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "s", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        break;
      case 'r':
        threshold=atof(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if((argc - optind) < 2)
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string train_f=argv[optind];
  std::string output=argv[optind+1];
  
  try
  {
    //tbb::task_scheduler_init init(2);
    string_table tbl;
    read_table(train_f.c_str(),tbl);
    if(!table_verify(tbl))
      return 1;
    
    //1. calculate averages
    volumes means;
    std::cout<<"Averaging..."<<std::endl;
    average_tables(tbl,means);
    std::cout<<"Done"<<std::endl;
    
    
    std::vector<gsl_double_matrix> cov(means.size());
    std::vector<gsl_double_matrix> pc(means.size());
    std::vector<gsl_double_matrix> ipc(means.size());
    std::vector<gsl_double_vector> w(means.size());
    std::vector<int> select(means.size());
    std::vector<double> weights(means.size());
    
    //2. PCA analysis for each modality
    for(int k=0;k<means.size();k++)
    {
      gsl_double_vector mod(tbl.size());
      mod=0;
      
      std::cout<<"Calculating covariance "<<k+1<<" out of "<<means.size()<<std::endl;
      calculate_covariance(tbl,means,cov[k],k,k+1,0,true);
      
      //calculating modulus of input vectors
      for(int i=0;i<mod.size();i++)
        mod.set(i,sqrt(cov[k].get(i,i))); 
        
      std::cout<<"Calculating PCA..."<<std::endl;
      select[k]=calculate_pca(cov[k],pc[k],w[k],threshold);
      std::cout<<"Selected "<<select[k]<<" components"<<std::endl;
      weights[k]=0.0;
      for(int i=0;i<select[k];i++) weights[k]+=fabs(w[k][i]);
      
      //2.5 calculate the individual weights
      _gsl_permutation perm(pc[k].rows());
      gsl_double_matrix lu;
      lu=pc[k];
      int sign;
      //calculating projection to the selected basis vectors
      for(int i=0;i<pc[k].rows();i++)
      {
        double norm=sqrt(fabs(w[k].get(i)));
        for(int j=0;j<pc[k].columns();j++)
        {
          lu.set(i,j,lu.get(i,j)/norm);
        }
      }
      gsl_linalg_LU_decomp(lu,perm,&sign);
      ipc[k].resize(pc[k].rows(),pc[k].columns());
      //normalize the vectors lengths
      gsl_linalg_LU_invert(lu,perm,ipc[k]);
      //now multiply by the sample's length
      for(int i=0;i<pc[k].rows();i++)
        for(int j=0;j<pc[k].columns();j++)
        {
          ipc[k].set(i,j,ipc[k].get(i,j)*mod.get(i));
        }
      std::cout<<std::endl;
    }
    
    //3. PCA analysis of the whole dataset
    gsl_double_matrix wcow;
    gsl_double_matrix wpc;
    gsl_double_vector ww;
    int wselect;
    
    //estimate weights of different modalities 
    //TODO: fix this to be more rasonable
    
    {
      double sum=0.0;
      for(int i=0;i<weights.size();i++)
        sum+=weights[i];

      for(int i=0;i<select.size();i++)
        weights[i]=weights[i]/sum;
    }
    
    calculate_covariance(ipc,select,weights,wcow);
    
    std::cout<<"Calculating overall PCA..."<<std::endl;
    wselect=calculate_pca(wcow,wpc,ww,threshold);
    std::cout<<"Weights of joint PCs:"<<std::endl;
    for(int i=0;i<ww.size();i++)
      std::cout<<ww[i]<<"\t";
    std::cout<<std::endl;
    std::cout<<"Selected "<<wselect<<" components"<<std::endl;
    
    //4. store results?
    std::cout<<"Writing basis vectors"<<std::endl;
    save_results(tbl,means,wpc,ww,output);
    std::cout<<"Done"<<std::endl;
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  
  return 0;
}
