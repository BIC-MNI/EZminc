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
#include <algorithm>

#include "pca_utils.h"
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <getopt.h>
#include <cmath>
#include <unistd.h>

#include <gsl_glue.h>
#include <gsl_gauss.h>

using namespace minc;

#if defined(_OPENMP)
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
#endif

void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" <training list> <output_prefix>  " << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask.mnc>"<<std::endl
      << "\t--buffers <n> default 1000 - number of buffers used for parallel processing"<<std::endl;
      //<< "\t--t-test output t-test results as well" << std::endl
}

bool run_fit(MNK_Gauss_Polinomial &mnk, 
             const std::vector< std::vector<double> > &design_matrix,
             const std::vector<double> &val,
             std::vector<double> &sol,
             double &RMS)
{
  mnk.clear();
  for(int i=0;i<design_matrix.size();i++)
    mnk.accumulate(design_matrix[i],val[i]);
  
  mnk.solve(sol);
  RMS=0.0;
  
  for(int i=0;i<design_matrix.size();i++)
  {
    double d=val[i]-mnk.fit(design_matrix[i],sol);
    RMS+=d*d;
  }
  RMS/=design_matrix.size()-1; //VF RMS is now == sd 
  RMS=sqrt(RMS);
  return true;
}


int main(int argc,char **argv)
{
  int clobber=0;
  int verbose=0;
  int normalize=0;
  int t_test=0;
  int parallel_buffers=1000;

  std::string mask_f;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"quiet",   no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"mask",    required_argument, 0, 'm'},
    {"buffers", required_argument, 0, 'b'},
    {"t-test",   no_argument, &t_test, 1},
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
      case 'b':
        parallel_buffers=atoi(optarg);
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
  std::string output_base=argv[optind+1];
  
  std::vector<std::string> output;
 
  try
  {

#ifdef _OPENMP
    std::cout<<"Using OpenMP, max number of threads="<<omp_get_max_threads()<<std::endl;
#endif
    
    
    string_table tbl;
    read_table_f(input_f.c_str(),tbl,1);
    if(!table_verify(tbl))
      return 1;
      
    //check if output exists
    
    int nout=tbl[0].size()-1;
    
    if(t_test)
      output.resize(nout*2+1);
     else
      output.resize(nout+1);
    
    for(int k=0;k<(nout+1);k++)
    {
      char tmp[1024];
      
      if(k<nout)
        sprintf(tmp,"%s_%d.mnc",output_base.c_str(),k);
      else
        sprintf(tmp,"%s_RMS.mnc",output_base.c_str()); 

      output[k]=tmp;
       
      if (!clobber && !access (output[k].c_str(), F_OK))
      {
        std::cerr << output[k].c_str () << " Exists!" << std::endl;
        return 1;
      }
      
      if(t_test&&k<nout)
      {
        sprintf(tmp,"%s_%d_t.mnc",output_base.c_str(),k);
        output[k*2+1]=tmp;
       
        if (!clobber && !access (output[k*2+1].c_str(), F_OK))
        {
          std::cerr << output[k*2+1].c_str () << " Exists!" << std::endl;
          return 1;
        }
      }
    }
    
    std::vector<minc_1_reader> files;
    std::vector<minc_input_iterator<float> > in;
    
    minc_1_reader mask;
    minc_input_iterator<unsigned char> mask_it;
      
    if(!mask_f.empty())
    {
      if(verbose)
        std::cout<<"Opening "<<mask_f.c_str()<<std::endl;
      mask.open(mask_f.c_str());
      mask.setup_read_byte();
      mask_it.attach(mask);
      mask_it.begin();
    }
    
    files.resize(tbl.size());
    in.resize(tbl.size());
      
    for(int i=0;i<tbl.size();i++)
    {
      if(verbose)
        std::cout<<"Opening "<<tbl[i][0].c_str()<<std::endl;
        
      files[i].open(tbl[i][0].c_str());
      
      if(i==0 && !mask_f.empty())
      {
        if(files[i].dim_no()!=mask.dim_no())
        {
          std::cerr<<"Input file "<< tbl[i][0].c_str() <<" should have same number of dimensions as mask!"<<std::endl;
          return 1;
        }
        bool good=true;
        for(int j=0;j<files[i].dim_no();j++)
          if(mask.dim(j).length!=files[i].dim(j).length)
            good=false;
        if(!good)
        {
          std::cerr<<"Input file "<< tbl[i][0].c_str() <<" should have same dimensions as mask!"<<std::endl;
          return 1;
        }
      }
      //check to make sure that all files are proper
      if(i>0)
      {
        if(files[0].dim_no()!=files[i].dim_no())
        {
          std::cerr<<"Input file "<< tbl[i][0].c_str() <<" should have same number of dimensions as first file!"<<std::endl;
          return 1;
        }
        bool good=true;
        for(int j=0;j<files[0].dim_no();j++)
          if(files[i].dim(j).length!=files[0].dim(j).length)
            good=false;
        if(!good)
        {
          std::cerr<<"Input file "<< tbl[i][0].c_str() <<" should have same dimensions as first file!"<<std::endl;
          return 1;
        }
      }
      files[i].setup_read_float();
      in[i].attach(files[i]);
      in[i].begin();
    }
    
    minc_info output_info=files[0].info();
    
    std::vector<minc_1_writer>                wrt(nout+1+(t_test?nout:0));
    std::vector<minc_output_iterator<float> > out(nout+1+(t_test?nout:0));
    
    if(verbose)
      std::cout<<"Writing to:";
    
    for(int k=0;k<(nout+1);k++)
    {
      if(verbose)
        std::cout<<output[k].c_str()<<" ";

      wrt[k].open(output[k].c_str(),output_info,2,NC_FLOAT);
      wrt[k].setup_write_float();

      out[k].attach(wrt[k]);
      out[k].begin();

      if(t_test&&k<nout)
      {
        if(verbose)
          std::cout<<output[k*2+1].c_str()<<" ";
      
        wrt[k*2+1].open(output[k*2+1].c_str(),output_info,2,NC_FLOAT);
        wrt[k*2+1].setup_write_float();
      
        out[k*2+1].attach(wrt[k]);
        out[k*2+1].begin();
      }
    }
    
    int len=files[0].ndim(1)*files[0].ndim(2)*files[0].ndim(3)*(files[0].ndim(0)?files[0].ndim(0):1)/parallel_buffers;
    
    int progress=0;
    int progress_steps=1;
    std::vector<int> report(100/progress_steps,0);
    std::vector< std::vector<double> > design_matrix(tbl.size());
    
    for(int i=0;i<tbl.size();i++)
    {
      design_matrix[i].resize(nout);
      for(int j=0;j<nout;j++)
      {
        design_matrix[i][j]=atof(tbl[i][j+1].c_str());
      }
    }
    

    std::vector<MNK_Gauss_Polinomial> mnk(parallel_buffers,MNK_Gauss_Polinomial(nout)); // linear approximation, we don't really use polinomial
    std::vector<std::vector<double> > val(parallel_buffers,std::vector<double>(tbl.size()));
    std::vector<std::vector<double> > sol(parallel_buffers,std::vector<double>(nout));
    std::vector<bool>   good(parallel_buffers);
    std::vector<double> RMS(parallel_buffers);
    
    while(!out[0].last())
    {
      int p_cnt=0;
      if(verbose)
      {
        progress++;
        int pct=(progress*100)/len;
        if(!report[pct/progress_steps] )
        {
          std::cout<<pct<<"% "<<std::flush;
          report[pct/progress_steps]=1;
        }
      }
      
      // fill up data
      for(p_cnt=0;p_cnt<parallel_buffers && !out[0].last(); p_cnt++)
      {
        if(!mask_f.empty())
        {
          if(mask_it.last()) break;
          good[p_cnt]=(mask_it.value()!=0);
          mask_it.next();
        } else {
          good[p_cnt]=true;
        }
        
        for(int i=0;i<tbl.size();i++)
        { 
          val[p_cnt][i]=in[i].value();
          in[i].next();
        }
      }
      
      //process data
      #pragma omp parallel for schedule(dynamic)
      for(int i=0;i<p_cnt;i++)
      {
        if(good[i] )
        {
          run_fit(mnk[i],design_matrix,val[i],sol[i],RMS[i]);
        } else {
          RMS[i]=0.0;
        }
      }
      
      // put it all back
      for(int i=0;i<p_cnt;i++)
      {
        for(int k=0;k<nout;k++)
        {
          out[k].value(sol[i][k]);
          out[k].next();
        }
        //last one is RMS
        out[nout].value(RMS[i]);
        out[nout].next();
      }
    }
    
    if(verbose) {
      std::cout<<"done!" << std::endl;
    }
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  
  return 0;
}
