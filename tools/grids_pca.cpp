/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2009 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include "minc_1_rw.h"
#include <iostream>
#include <fstream>

#include "pca_utils.h"
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <getopt.h>
#include <cmath>
#include <gsl_glue.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

using namespace minc;

void save_results(const string_table& tbl, 
                  const grids &mean, 
                  const minc_byte_volume& mask,
                  gsl_double_matrix &pc, 
                  gsl_double_vector &w, 
                  const std::string &output,
                  bool normalize,
                  double dv)
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
      save_simple_volume(wrt,mean[k]);
    }
    
    minc_grid_volume tmp,tmp_read;
    
    for(int j=0;j<pc.size(0);j++)
    {
      char file2[1024];
      sprintf(file2,"%s_%d_%03d.mnc",output.c_str(),k,j+1); // 0 - mean
      tmp.resize(mean[k].size());
      out_tbl[j+1][k+1]=minc::_basename(file2);
      char tmp2[256];
      sprintf(tmp2,"%g",w[j]);
      out_tbl[j+1][0]=tmp2;
      tmp=IDX<float>(0.0,0.0,0.0);
      for(int i=0;i<pc.size(1);i++)
      {
        minc_1_reader rdr2;
        rdr2.open(tbl[i][k].c_str());
        load_simple_volume(rdr2,tmp_read);
        tmp_read-=mean[k];
        tmp.weighted_add(tmp_read,pc.get(i,j));
      }
      minc_1_writer wrt2; 
      wrt2.open(file2,rdr.info(),2,NC_FLOAT);
      
      if(normalize)
      {
        float len=sqrt(sum2(tmp,mask)*dv);
        tmp/=IDX(len,len,len);
      }
      save_simple_volume(wrt2,tmp);
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


void save_results_cached(const string_table& tbl, 
                  const grids &mean,
                  const minc_byte_volume& mask,
                  gsl_double_matrix &pc, 
                  gsl_double_vector &w, 
                  const std::string &output,
                  bool normalize, double dv)
{

  string_table out_tbl(pc.size(0)+1);
  
  for(int j=0;j<(pc.size(0)+1);j++)
  {
    out_tbl[j].resize(tbl[0].size()+1);
  }
  
  std::vector<grids> inputs;
  inputs.resize(tbl.size());
  for(int j=0;j<tbl.size();j++)
  {
    inputs[j].resize(tbl[0].size());
    for(int i=0;i<tbl[0].size();i++)
    {
      minc_1_reader rdr;
      rdr.open(tbl[j][i].c_str());
      load_simple_volume(rdr,inputs[j][i]);
      inputs[j][i]-=mean[i];
    }
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
      save_simple_volume(wrt,mean[k]);
    }
    
    minc_grid_volume tmp,tmp_read;
    
    for(int j=0;j<pc.size(0);j++)
    {
      char file2[1024];
      sprintf(file2,"%s_%d_%03d.mnc",output.c_str(),k,j+1); // 0 - mean
      tmp.resize(mean[k].size());
      out_tbl[j+1][k+1]=minc::_basename(file2);
      char tmp2[256];
      sprintf(tmp2,"%g",w[j]);
      out_tbl[j+1][0]=tmp2;
      tmp=IDX<float>(0.0,0.0,0.0);
      for(int i=0;i<pc.size(1);i++)
      {
        tmp.weighted_add(inputs[i][k],pc.get(i,j));
      }
      
      if(normalize)
      {
        float len=sqrt(sum2(tmp,mask)*dv);
        tmp/=IDX(len,len,len);
      }
      minc_1_writer wrt2;
      wrt2.open(file2,rdr.info(),2,NC_FLOAT);
      save_simple_volume(wrt2,tmp);
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
      << "\t--mask <mask.mnc> "                 << std::endl
      << "\t--threshold <f> threshold for the principal components, default 0.98" <<std::endl
      << "\t--nomean don't substract mean" <<std::endl
      << "\t--cache  - cache all files in memory"<<std::endl
      << "\t--normalize - normalize the lengths of principal components to 1.0"<<std::endl
      ;
}

int main(int argc,char **argv)
{
  int clobber=0;
  int verbose=0;
  double threshold=0.98;
  std::string mask_f;
  int cache=0;
  int nm=0;
  int normalize=0;
  static struct option long_options[] =
  {
    {"verbose",   no_argument,       &verbose,  1},
    {"quiet",     no_argument,       &verbose,  0},
    {"normalize",   no_argument,       &normalize,  1},
    {"nonormalize", no_argument,       &normalize,  0},
    {"cache",     no_argument,       &verbose,  1},
    {"clobber",   no_argument,       &clobber,  1},
    {"threshold", required_argument, 0,       'r'},
    {"mask"     , required_argument, 0,       'm'},
    {"nomean"   , no_argument,       &nm,       1},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "rm:", long_options, &option_index);

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

  if((argc - optind) < 2)
  {
    show_usage(argv[0]);
    return 1;
  }
  
  std::string train_f=argv[optind];
  std::string output=argv[optind+1];
  std::string output_rotation=output+"_rotation.csv";
  
  try
  {
    double dv=1;
    string_table tbl;
    read_table(train_f.c_str(),tbl);
    if(!table_verify(tbl))
      return 1;
    
    dv=calc_dv(tbl[0][0].c_str());
    
    //1. calculate averages
    grids means;
    if(nm)
    {
      means.resize(tbl[0].size());
      for(int i=0;i<tbl[0].size();i++)
      {
        minc_1_reader rdr;
        rdr.open(tbl[0][i].c_str(),false);
        means[i].resize(rdr.ndim(1),rdr.ndim(2),rdr.ndim(3));
        means[i]=IDX<float>(0,0,0); //0 mean
      }
    } else {
      std::cout<<"Averaging..."<<std::endl;
      average_tables(tbl,means);
      std::cout<<"Done"<<std::endl;
    }
    minc_byte_volume mask;
    if(!mask_f.empty())
    {
      minc_1_reader rdr;
      rdr.open(mask_f.c_str());
      load_simple_volume(rdr,mask);
    } else {
      mask.resize(means[0].size());
      mask=1; //all is allowed
    }
    
    gsl_double_matrix cov;
    gsl_double_matrix pc;
    gsl_double_vector w;
    int select;
    
    std::cout<<"Calculating covariance "<<std::endl;
    if(cache)
      calculate_covariance_cached(tbl,means,mask,cov,0,means.size(),0,dv);
    else
      calculate_covariance(tbl,means,mask,cov,0,means.size(),0,dv);

    std::cout<<"Calculating PCA..."<<std::endl;
    select=calculate_pca(cov,pc,w,threshold);
    std::cout<<"Selected "<<select<<" components out of "<<tbl.size()<<std::endl;

    std::ofstream out_r(output_rotation.c_str());

    for(int j=0;j<pc.size(1);j++)
    {
      for(int i=0;i<pc.size(0);i++)
      {
        out_r<<pc.get(i,j)<<",";
      }
      out_r<<std::endl;
    }

    //4. store results
    if(verbose) 
      std::cout<<"Writing basis vectors"<<std::endl;

    if(cache)
      save_results_cached(tbl,means,mask,pc,w,output,normalize,dv);
    else
      save_results(tbl,means,mask,pc,w,output,normalize,dv);

    std::cout<<"Done"<<std::endl;

  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  return 0;
}
