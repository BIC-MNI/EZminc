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
#include <gsl_glue.h>
#include <gsl_gauss.h>
#include <gsl/gsl_blas.h>
#include <strtok.h>
#include <unistd.h>

using namespace minc;

typedef std::vector<gsl_double_vector> gsl_vectors;
typedef std::vector<gsl_double_matrix> gsl_matrixes;


//! print mu & sigma
void print_mu_sigma( gsl_vectors&      mu, gsl_matrixes&     sigma)
{
  int number_of_classes=mu.size();
  int number_of_inputs=mu[0].size();
  
  std::cout<<"Means:"<<std::endl;
  for(int i=0;i<number_of_classes;i++)
  {
    if(i)
      std::cout<<"\t ";
    else
      std::cout<<"\t[";
    for(int k=0;k<number_of_inputs;k++)
    {
      if(k) std::cout<<",";
      std::cout<<mu[i].get(k);
    }
    if(i!=(number_of_classes-1))
      std::cout<<";"<<std::endl;
    else
      std::cout<<"]"<<std::endl;
  }
      
  std::cout<<"Sigmas:"<<std::endl;
      
  for(int i=0;i<number_of_classes;i++)
  {
    for(int k=0;k<number_of_inputs;k++)
    {
      if(k)
        std::cout<<"\t ";
      else 
        std::cout<<"\t[ ";
      for(int j=0;j<number_of_inputs;j++)
      {
        if(j) std::cout<<",";
        std::cout<<sigma[i].get(k,j);
      }
      if(k!=(number_of_inputs-1))
        std::cout<<";"<<std::endl;
      else
        std::cout<<"]"<<std::endl;
    }
  }
}

//! estimate sample mu & sigma using discrete classes
void estimate_mu_and_sigma(volumes & input,
                           minc_byte_volume& mask,
                           minc_byte_volume& cls,
                           gsl_vectors&      mu,
                           gsl_matrixes&     sigma,
                           double &small)
{
  int number_of_inputs=input.size();
  int number_of_classes=mu.size();
  int i,j,k;
  gsl_vectors  _means(number_of_classes,gsl_double_vector(number_of_inputs,0));
  gsl_matrixes _sigma(number_of_classes,gsl_double_matrix(number_of_inputs,number_of_inputs,0.0));
  std::vector<int>    _counts(number_of_classes,0);
  
  std::vector<double> _total_means(number_of_inputs,0.0);
  std::vector<double> _total_sigma(number_of_inputs,0.0);
  
  int _count=0;
  
  //1 calculate means
  for(j=0;j<cls.c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) continue; //only use ROI
    _count++;
    
    int _c=cls.c_buf()[j];
    if(!_c || _c>number_of_classes) continue; //only use classified voxel
    _c--;
    _counts[_c]++;
    
    for(i=0;i<number_of_inputs;i++)
    {
      
      double v=input[i].c_buf()[j];
      _total_means[i]+=v;
      _total_sigma[i]+=v*v;
      
      _means[_c].set(i,_means[_c].get(i)+v);
    }
  }
  
  if(!_count)
    REPORT_ERROR("No voxels defined in ROI!");
  
  for(k=0;k<number_of_classes;k++)
  {
    
    for(i=0;i<number_of_inputs;i++)
      if(_counts[k]>0)
        _means[k].set(i,_means[k].get(i)/_counts[k]);
      else
        _means[k].set(i,mu[k].get(i));
  }
  
  for(i=0;i<number_of_inputs;i++)
  {
    _total_means[i]/=_count/number_of_inputs;
    _total_sigma[i]/=_count/number_of_inputs;
    _total_sigma[i]-=_total_means[i]*_total_means[i];
  }
  
  //2 estimate sigmas
  for(j=0;j<cls.c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) continue; //only use ROI
    int _c=cls.c_buf()[j];
    if(!_c || _c>number_of_classes) continue; //only use classified voxel
    _c--;
    
    for(i=0;i<number_of_inputs;i++)
    {
      double v1=input[i].c_buf()[j]-_means[_c].get(i);
      
      for(int t=0;t<=i;t++)
      {
        double v2=input[t].c_buf()[j]-_means[_c].get(t);
        _sigma[_c].set(i,t,_sigma[_c].get(i,t)+v1*v2);
      }
    }
  }
  
  for(k=0;k<number_of_classes;k++)
  {
    for(i=0;i<number_of_inputs;i++)
    {
      for(int t=0;t<=i;t++)
      {
        if(_counts[k]>0)
          _sigma[k].set(i,t,_sigma[k].get(i,t)/_counts[k]);
        else 
        {
          _sigma[k].set(i,t,sigma[k].get(i,t));
          if(i==t && _sigma[k].get(i,t)==0.0) //assume this is not set
          {
            _sigma[k].set(i,i,_total_sigma[i]/(number_of_classes*number_of_classes));
          }
        }
        
        sigma[k].set(i,t,_sigma[k].get(i,t));
        
        if(i!=t)
          sigma[k].set(t,i,_sigma[k].get(i,t));
      }
      mu[k].set(i,_means[k].get(i));
    }
  }
  
  // calculating tolerances:
  small=0;
  for(i=0;i<number_of_inputs;i++)
    small+=_total_sigma[i];
  
  small/=number_of_inputs*number_of_classes;
}


//! estimate sample mu & sigma using given PDFs
void estimate_mu_and_sigma(volumes & input,
                           minc_byte_volume& mask,
                           volumes& p,
                           gsl_vectors&      mu,
                           gsl_matrixes&     sigma,
                           std::vector<double>& priors)
{
  int number_of_inputs=input.size();
  int number_of_classes=mu.size();
  int i,j,k;
  
  gsl_vectors  _means(number_of_classes,gsl_double_vector(number_of_inputs,0.0));
  std::vector<double>  _sums(number_of_classes,0.0);
  
  gsl_matrixes _sigma(number_of_classes,gsl_double_matrix(number_of_inputs,number_of_inputs,0.0));
  int _count=0;
  
  //1 calculate means (weighted sum)
  for(j=0;j<mask.c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) continue; //only use ROI
    _count++;
    
    for(k=0;k<number_of_classes;k++)
    {
      _sums[k]+=p[k].c_buf()[j];
      
      for(i=0;i<number_of_inputs;i++)
      {
        _means[k].set(i,_means[k].get(i)+input[i].c_buf()[j]*p[k].c_buf()[j]);
      }
    }
  }
  
  for(k=0;k<number_of_classes;k++)
  {
    for(i=0;i<number_of_inputs;i++)
      if(_sums[k]>0)
        _means[k].set(i,_means[k].get(i)/_sums[k]);
      else //no mean for this class?
        _means[k].set(i,mu[k].get(i));
  }
  
  //2 estimate sigmas
  for(j=0;j<mask.c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) continue; //only use ROI
    
    for(k=0;k<number_of_classes;k++)
    {
      if(_sums[k]==0.0) 
        continue;
      
      for(i=0;i<number_of_inputs;i++)
      {
        double v1=input[i].c_buf()[j]-_means[k].get(i);
        
        for(int t=0;t<=i;t++)
        {
          double v2=input[t].c_buf()[j]-_means[k].get(t);
          
          _sigma[k].set(i,t,_sigma[k].get(i,t)+v1*v2*p[k].c_buf()[j]);
        }
      }
    }
  }
  
  for(k=0;k<number_of_classes;k++)
  {
    if(_sums[k]>0)
    {
      for(i=0;i<number_of_inputs;i++)
      {
        for(int t=0;t<=i;t++)
        {
          sigma[k].set(i,t,_sigma[k].get(i,t)/_sums[k]);
          
          if(i!=t)
            sigma[k].set(t,i,sigma[k].get(i,t));
        }
        mu[k].set(i,_means[k].get(i));
      }
      priors[k]=_sums[k]/_count;
    }
  }
}


//! calculate PDFs
void calculate_probabilities(volumes & input,
                            minc_byte_volume& mask,
                            gsl_vectors&      mu,
                            gsl_matrixes&     sigma,
                            double small,
                            volumes & p
                            )
{
  int number_of_inputs=input.size();
  int number_of_classes=mu.size();
  int i,j,k;
  
  double vsmall=small/1000;
  
  if(p.size()!=  number_of_classes)
  {
    p.resize(number_of_classes,input[0]);
  }
  //calculate invert sigma
  gsl_matrixes U(sigma);
  gsl_matrixes V(number_of_classes, gsl_double_matrix(number_of_inputs,number_of_inputs,0));
  gsl_vectors  S(number_of_classes, gsl_double_vector(number_of_inputs,0));
  gsl_double_vector svd_work(number_of_inputs,0);
  
  std::vector<double> detSigma(number_of_classes,1.0);
  std::vector<double> dims(number_of_classes,0);
  
  for(k=0;k<number_of_classes;k++)
  {
    gsl_linalg_SV_decomp_jacobi(U[k],V[k],S[k]);
    for(i=0;i<number_of_inputs;i++)
    {
      if(fabs(S[k].get(i))<vsmall) 
        S[k].set(i,0.0);
      else
      {
        detSigma[k]*=fabs(S[k].get(i));
        dims[k]++;
      }
      
    }
    detSigma[k]=sqrt(detSigma[k]);
  }
  
  gsl_double_vector _v(number_of_inputs),_x(number_of_inputs);
  
  //calculate all the classes
  for(k=0;k<number_of_classes;k++)
  {
    double coeff=1.0/(pow(2*M_PI,dims[k]/2.0)*detSigma[k]);
    for(j=0;j<input[0].c_buf_size();j++)
    {
      if(mask.c_buf()[j]) 
      {
        for(int i=0;i<number_of_inputs;i++)
        {
          _v.set(i,input[i].c_buf()[j]-mu[k].get(i));
        }
        gsl_linalg_SV_solve(U[k],V[k],S[k],_v,_x);
        double prod;
        gsl_blas_ddot(_x,_v,&prod);
        p[k].c_buf()[j]=coeff*exp(-prod/2.0);
      } else {
        p[k].c_buf()[j]=0.0;
      }
    }
  }
}

//! estimate sample mu & sigma using given PDFs
void estimate_mu_and_sigma_interpolate(volumes & input,
                           minc_byte_volume& mask,
                           volumes& p,
                           gsl_vectors&      mu,
                           gsl_matrixes&     sigma,
                           std::vector<double>& priors,
                           int steps=1)
{
  int number_of_inputs=input.size();
  int number_of_classes=mu.size();
  int i,j,k;
  
  gsl_vectors  _means(number_of_classes,gsl_double_vector(number_of_inputs,0.0));
  std::vector<double>  _sums(number_of_classes,0.0);
  
  gsl_matrixes _sigma(number_of_classes,gsl_double_matrix(number_of_inputs,number_of_inputs,0.0));
  int _count=0;
  
  //1 calculate means (weighted sum) (trilinear interpolation shouldn's change means)
  for(j=0;j<mask.c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) continue; //only use ROI
    _count++;
    
    for(k=0;k<number_of_classes;k++)
    {
      _sums[k]+=p[k].c_buf()[j];
      
      for(i=0;i<number_of_inputs;i++)
      {
        _means[k].set(i,_means[k].get(i)+input[i].c_buf()[j]*p[k].c_buf()[j]);
      }
    }
  }
  
  for(k=0;k<number_of_classes;k++)
  {
    for(i=0;i<number_of_inputs;i++)
      if(_sums[k]>0)
        _means[k].set(i,_means[k].get(i)/_sums[k]);
      else //no mean for this class?
        _means[k].set(i,mu[k].get(i));
  }
  
  //2 estimate sigmas (interpolate)
  //for(j=0;j<mask.c_buf_size();j++)
  for(int z=0;z<mask.dim(2);z++)
    for(int y=0;y<mask.dim(1);y++)
      for(int x=0;x<mask.dim(0);x++)
  
  {
    if(!mask.get(x,y,z)) continue; //only use ROI
    
    for(k=0;k<number_of_classes;k++)
    {
      if(_sums[k]==0.0) 
        continue;
      
      for(i=0;i<number_of_inputs;i++)
      {
        for(int t=0;t<=i;t++)
        {
          for(int zz=-steps;zz<=steps;zz++)
            for(int yy=-steps;yy<=steps;yy++)
              for(int xx=-steps;xx<=steps;xx++)
              {
                double v1=input[i].interpolate(x+(float)xx/(2.0*steps),y+(float)yy/(2.0*steps),z+(float)zz/(2.0*steps))
                   -_means[k].get(i);
                
                double v2=input[t].interpolate(x+(float)xx/(2.0*steps),y+(float)yy/(2.0*steps),z+(float)zz/(2.0*steps))
                    -_means[k].get(t);
                
                _sigma[k].set(i,t,_sigma[k].get(i,t)+v1*v2*p[k].get(x,y,z));
              }
        }
      }
    }
  }
  
  for(k=0;k<number_of_classes;k++)
  {
    if(_sums[k]>0)
    {
      for(i=0;i<number_of_inputs;i++)
      {
        for(int t=0;t<=i;t++)
        {
          sigma[k].set(i,t,_sigma[k].get(i,t)/(_sums[k]*(steps*2+1)*(steps*2+1)*(steps*2+1)));
          
          if(i!=t)
            sigma[k].set(t,i,sigma[k].get(i,t));
        }
        mu[k].set(i,_means[k].get(i));
      }
      priors[k]=_sums[k]/_count;
    }
  }
}


//! calculate PDFs
void calculate_probabilities_interpolate(volumes & input,
                            minc_byte_volume& mask,
                            gsl_vectors&      mu,
                            gsl_matrixes&     sigma,
                            double small,
                            volumes & p,
                            int steps=1 )
{
  int number_of_inputs=input.size();
  int number_of_classes=mu.size();
  int i,j,k;
  
  double vsmall=small/1000;
  
  if(p.size()!=  number_of_classes)
  {
    p.resize(number_of_classes,input[0]);
  }
  //calculate invert sigma
  gsl_matrixes U(sigma);
  gsl_matrixes V(number_of_classes, gsl_double_matrix(number_of_inputs,number_of_inputs,0));
  gsl_vectors  S(number_of_classes, gsl_double_vector(number_of_inputs,0));
  gsl_double_vector svd_work(number_of_inputs,0);
  
  std::vector<double> detSigma(number_of_classes,1.0);
  std::vector<double> dims(number_of_classes,0);
  
  for(k=0;k<number_of_classes;k++)
  {
    gsl_linalg_SV_decomp_jacobi(U[k],V[k],S[k]);
    for(i=0;i<number_of_inputs;i++)
    {
      if(fabs(S[k].get(i))<vsmall) 
        S[k].set(i,0.0);
      else
      {
        detSigma[k]*=fabs(S[k].get(i));
        dims[k]++;
      }
      
    }
    detSigma[k]=sqrt(detSigma[k]);
  }
  
  gsl_double_vector _v(number_of_inputs),_x(number_of_inputs);
  
  //calculate all the classes
  for(k=0;k<number_of_classes;k++)
  {
    double coeff=1.0/(pow(2*M_PI,dims[k]/2.0)*detSigma[k]);
    //for(j=0;j<input[0].c_buf_size();j++)
    for(int z=0;z<mask.dim(2);z++)
      for(int y=0;y<mask.dim(1);y++)
        for(int x=0;x<mask.dim(0);x++)
    
    {
      p[k].set(x,y,z,0);
      if(mask.get(x,y,z)) 
      {
        //calculate average P ?
        //calculating max p //TODO replace with log-exp average
        double mp=0;
        for(int zz=-steps;zz<=steps;zz++)
          for(int yy=-steps;yy<=steps;yy++)
            for(int xx=-steps;xx<=steps;xx++)
            {
              for(int i=0;i<number_of_inputs;i++)
              {
                _v.set(i,input[i].interpolate(x+(float)xx/(2.0*steps),y+(float)yy/(2.0*steps),z+(float)zz/(2.0*steps))-mu[k].get(i));
              }
              gsl_linalg_SV_solve(U[k],V[k],S[k],_v,_x);
              double prod;
              gsl_blas_ddot(_x,_v,&prod);
              double _p=coeff*exp(-prod/2.0);
              if(_p>mp) mp=_p;
            }
        p[k].set(x,y,z,mp);
      } 
    }
  }
}


//! apply uniform apriory information
void apply_priors(volumes & p,minc_byte_volume& mask,std::vector<double>& priors)
{
  int number_of_classes=p.size();
  
  for(int j=0;j<p[0].c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) continue;
    
    for(int k=0;k<number_of_classes;k++)
    {
      p[k].c_buf()[j]*=priors[k];
    }
  }
}

//! apply uniform apriory information
void apply_priors(volumes & p,minc_byte_volume& mask,volumes& priors)
{
  int number_of_classes=p.size();
  
  for(int j=0;j<p[0].c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) continue;
    
    for(int k=0;k<number_of_classes;k++)
    {
      double v=priors[k].c_buf()[j];
      if(v<0) v=0;
      p[k].c_buf()[j]*=v;
    }
  }
}

//!normalize probabilities so that total is equal to 1 at each voxel
void normalize_probabilities(volumes & p,minc_byte_volume& mask)
{
  int number_of_classes=p.size();
  
  int cnt=0;
  
  for(int j=0;j<p[0].c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) 
    {
      for(int k=0;k<number_of_classes;k++)
        p[k].c_buf()[j]=0.0;
      continue;
    }
    cnt++;
    
    double sum=0.0;
    for(int k=0;k<number_of_classes;k++)
    {
      sum+=p[k].c_buf()[j];
    }
    
    if(sum>0) //we have something here
      for(int k=0;k<number_of_classes;k++)
        p[k].c_buf()[j]/=sum;
    else
      for(int k=0;k<number_of_classes;k++)
        p[k].c_buf()[j]=1.0/number_of_classes;
  }
}

void apply_hard_classify(volumes & p,minc_byte_volume& mask,minc_byte_volume& cls)
{
  int number_of_classes=p.size();
  
  if(cls.size()!=mask.size())
    cls.resize(mask.size());
  
  for(int j=0;j<p[0].c_buf_size();j++)
  {
    if(!mask.c_buf()[j]) 
    {
      cls.c_buf()[j]=0;
      continue;
    }
    double max_p=0.0;
    int c=0;
    for(int k=0;k<number_of_classes;k++)
    {
      if(p[k].c_buf()[j]>max_p) 
      {
        max_p=p[k].c_buf()[j];
        c=k;
      }
    }
    c++;
    cls.c_buf()[j]=c;
  }
}

void show_usage(const char *name)
{
  std::cerr 
      << "Usage: "<<name<<" <input_1> .. <input_n> <output>  " << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--mask <mask.mnc> ROI mask"<<std::endl
      << "\t--train <train.mnc> training file"<<std::endl
      << "\t--means <c1,c2...,cn> precomputed class means"<<std::endl
      << "\t--sigma <s1,.....,sn> precomputed class sigmas"<<std::endl
      << "\t--classes <n> number of classes"<<std::endl
      << "\t--prob p1,p2... apriori probabilites of each class"<<std::endl
      << "\t--iter <n> maximum number of iterations"<<std::endl
      << "\t--priors <p1.mnc>,<p2.mnc>,...,<pn.mnc>"<<std::endl
      << "\t--save_prob <base> save probability maps"<<std::endl
      << "\t--debug print more info"<<std::endl
      << "\t--interpolate use intensity interolation"<<std::endl;
}

int main(int argc,char **argv)
{
  int clobber=0;
  int verbose=0;
  int normalize=0;
  int debug=0;
  std::string prob_f;
  std::string mask_f;
  std::string train_f;
  int number_of_classes=2;
  std::vector<double> _means;
  std::vector<double> _sigmas;
  std::vector<double> _probs;
  std::vector<std::string> _priors_f;
  int maxiter=10;
  int interpolate=0;
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"debug",   no_argument, &debug, 1},
    {"quiet",   no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"interpolate", no_argument, &interpolate, 1},
    {"mask",    required_argument,    0, 'm'},
    {"train",    required_argument,   0, 't'},
    {"means",    required_argument,   0, 'M'},
    {"sigma",    required_argument,   0, 's'},
    {"classes",  required_argument, 0, 'c'},
    {"prob",    required_argument, 0, 'p'},
    {"iter",    required_argument, 0, 'i'},
    {"priors",  required_argument, 0, 'r'},
    {"save_prob",  required_argument, 0, 'S'},
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
      case 't':
        train_f=optarg;
        break;
      case 'M':
        string_to_vector(optarg,_means);
        break;
      case 's':
        string_to_vector(optarg,_sigmas);
        break;
      case 'p':
        string_to_vector(optarg,_probs);
        break;
      case 'r':
        stringtok(_priors_f,optarg,", ");
        break;
      case 'c':
        number_of_classes=atoi(optarg);
        break;
      case 'i':
        maxiter=atoi(optarg);
        break;
      case 'S':
        prob_f=optarg;
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

  std::string output_f=argv[argc-1]; //last is output
  argc--;

  strings _input_f;
  for(int i=optind;i<argc;i++)
    _input_f.push_back(std::string(argv[i]));
  
  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr<<"File "<<output_f.c_str()<<" exists, use --clobber"<<std::endl;
    return 1;
  }
  
  strings _outputs_f;
  
  if(!prob_f.empty())
  {
    for(int i=0;i<number_of_classes;i++)
    {
      char tmp[1024];
      sprintf(tmp,"%s_%d.mnc",prob_f.c_str(),i+1);
      _outputs_f.push_back(tmp);

      if (!clobber && !access (tmp, F_OK))
      {
        std::cerr<<"File "<<tmp<<" exists, use --clobber"<<std::endl;
        return 1;
      }
    }
  }
      
  int number_of_inputs=_input_f.size();

  try
  {
    if(!_means.empty() && _means.size()!=number_of_classes*number_of_inputs)
      REPORT_ERROR("Incorrect number of means");
    
    if(!_sigmas.empty() && _sigmas.size()!=number_of_classes*number_of_inputs)
      REPORT_ERROR("Incorrect number of sigmas");
      
    if(!_probs.empty() && _probs.size()!=number_of_classes)
      REPORT_ERROR("Incorrect number of apriori probabilites");
    
    if(!_priors_f.empty() && _priors_f.size()!=number_of_classes)
      REPORT_ERROR("Incorrect number of priors");
    
    if(train_f.empty() && _means.empty()) 
      REPORT_ERROR("Specify either training segmentation or means");
    
    //1  load all the volumes
    volumes input,prob,priors;
    if(verbose) std::cout<<"loading input files"<<std::endl;
    
    load_volumes(_input_f,input,0); 
    check_volumes(input);
    
    if(!_priors_f.empty())
    {
      load_volumes(_priors_f,priors,0);
      check_volumes(priors);
      if(priors[0].size()!=input[0].size())
        REPORT_ERROR("Priors volume size mismatch!");
    }
    std::cout<<priors.size()<<std::endl;
    
    minc_byte_volume mask,cls;
    
    if(!mask_f.empty())
    {
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<mask_f.c_str()<<std::endl;
      rdr.open(mask_f.c_str());
      load_simple_volume(rdr,mask);
      if(mask.size()!=input[0].size())
        REPORT_ERROR("Mask size mismatch");
    } else {
      mask.resize(input[0].size());
      mask=1; //all volume is selected
    }
    
    if(!train_f.empty())
    {
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading training data:"<<train_f.c_str()<<std::endl;
      rdr.open(train_f.c_str());
      load_simple_volume(rdr,cls);
      if(cls.size()!=input[0].size())
        REPORT_ERROR("Training  size mismatch");
    } else {
      cls.resize(input[0].size());
      cls=0; //none voxel is selected
    }
    gsl_vectors     mu(number_of_classes);
    gsl_matrixes sigma(number_of_classes);
    
    for(int i=0;i<number_of_classes;i++)
    {
      mu[i].resize(input.size());
      mu[i].set_all(0);
      sigma[i].resize(input.size(),input.size());
      sigma[i].set_all(0);
      
      if(!_means.empty())
      {
        for(int j=0;j<input.size();j++)
          mu[i].set(j,_means[j+i*input.size()]);
      }
      
      if(!_sigmas.empty())
      {
        for(int j=0;j<input.size();j++)
          sigma[i].set(j,j,_sigmas[j+i*input.size()]);
      }
    }
    
    if(_probs.empty())
      _probs.resize(number_of_classes,1.0/number_of_classes);
   
    // make initial estimate 
    double small=1e-5;
    if(!train_f.empty())
    {
      if(verbose) std::cout<<"Estimating initial parameters"<<std::endl;
      estimate_mu_and_sigma(input,mask,cls,mu,sigma,small);
    }
    
    
    for(int iter=0;iter<maxiter;iter++)
    {
      if(debug)
      {
        print_mu_sigma(mu,sigma);
      }
      
      if(verbose)
      {
        std::cout<<"Prob: ";
        for(int k=0;k<number_of_classes;k++)
          std::cout<<_probs[k]<<" ";
        std::cout<<std::endl;
          
        std::cout<<"Calculating probabilities (E) "<<iter<<std::endl;
      }
      
      if(interpolate) 
        calculate_probabilities_interpolate(input,mask,mu,sigma,small,prob);
      else
        calculate_probabilities(input,mask,mu,sigma,small,prob);
        
      if(debug)
      {
        print_mu_sigma(mu,sigma);
      }
      
      if(priors.empty())
        apply_priors(prob, mask,_probs);
      else
        apply_priors(prob, mask,priors);
      
      normalize_probabilities(prob,mask);
      if(verbose) std::cout<<"Estimating new parameters (M) "<<std::endl;
      
      if(iter==(maxiter-1)) break; //no need to estimate parameters on last loop
      if(interpolate)
        estimate_mu_and_sigma_interpolate(input,mask,prob,mu,sigma,_probs);
      else
        estimate_mu_and_sigma(input,mask,prob,mu,sigma,_probs);
    }
    
    apply_hard_classify(prob,mask,cls);
    save_volume(cls,output_f.c_str(),_input_f[0].c_str());
    
    if(!prob_f.empty())
    {
      if(verbose) std::cout<<"Saving probabilities"<<std::endl;
      save_volumes(_outputs_f,prob,_input_f[0].c_str());
    }
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
}
