/* ----------------------------- MNI Header -----------------------------------
@NAME       :  fit_harmonics_grids
@DESCRIPTION:  spherical harminic LSQ approximation programm
@COPYRIGHT  :
              Copyright 2006 Vladimir Fonov, McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>

#include "gsl_glue.h"
#include "gsl_gauss.h"

#include <itkBSplineInterpolateImageFunction.h>
#include <unistd.h>
#include <getopt.h>

#include "sphericalHarmonicsTransform.h"
#include <gsl/gsl_rng.h>
#include <time.h>
#include <sys/time.h>

#include "minc_wrappers.h"


#ifdef HAVE_MINC4ITK
#include <itkMincGeneralTransform.h>
typedef minc::XfmTransform<double,3> XFMTransformType;
#else
#include <itkMINCTransformAdapter.h>
typedef itk::MINCTransformAdapter<double,3,3> XFMTransformType;
#endif


using namespace std;
using namespace minc;

void print_coeff (std::ostream & out, const std::vector < double >&coeff2)
{
  out.precision(40);
  int order = coeff2.size ();
  for (int i = 0; i < order; i++) {
    out << coeff2[i] << " ";
    if(!out.good())
      REPORT_ERROR("Can't write to file");
  }
  out << std::endl;
}

void save_coeff (std::ostream & out, const std::vector < double >&coeff2)
{
  int order = coeff2.size ();
  out.precision(40);
  for (int i = 0; i < order; i++) {
    out << coeff2[i] << endl;
    if(!out.good())
      REPORT_ERROR("Can't write to file");
  }
}

class Tag_fit
{
protected:
  double _last_distance;
  static const double _distance_epsilon;
public:
  typedef std::vector <double> fitting_coeff;
  typedef std::vector <fitting_coeff> fittings;
  typedef std::vector <int>    Index;
  typedef std::vector <double> Distances;

  int       order;
  double    keep;
  bool      verbose;
  double    max_dev;
  int       _max_iterations;
  double    sd,max_distance;
  fittings  coeff;
  bool      limit_linear;
  int       skip_voxels;
  bool      use_svd;
  double    scaling;

  //fittings  basis_x,basis_y,basis_z;

  Index     index;
  Distances distances;
  double    lambda;

  std::vector<minc::mask3d::Pointer>        _mask;
  std::vector<XFMTransformType::Pointer>    _lin_xfm;
  std::vector<XFMTransformType::Pointer>    _nl_xfm;
  
  size_t  _sample_count;

  basis_functions_x fun_x; // here we have the same basis for X,Y,Z

  Tag_fit(int order_, double keep_,
          double max_dev_=5.0,
          int max_iter=100,bool verbose_=false,
          bool ll=false,int skip_voxels_=0):
    order(order_),
    keep(keep_),
    verbose(verbose_),
    max_dev(max_dev_),
    _max_iterations(max_iter),
    sd(0),
    max_distance(0),
    coeff(3),
    limit_linear(ll),
    skip_voxels(skip_voxels_),
    lambda(0.5),_sample_count(0),use_svd(true),
    scaling(200.0)
  {
  }

  void fit_coeff(bool condition=false) {
    coeff.resize(3);
    coeff[0].resize(order);
    coeff[1].resize(order);
    coeff[2].resize(order);

    basis_vector bas_x1(order);
    basis_vector bas_x2(order);
    
    fun_x.set_scaling(scaling);

#ifdef SEPARATE    
    minc::LSQ_Gauss<double> pol[3];
    pol[0].resize(order-3);
    pol[1].resize(order-3);
    pol[2].resize(order-3);
    basis_vector bas_diff(order-3);
#else
    minc::LSQ_Gauss<double> pol((order-3)*3);
    basis_vector bas_diff((order-3)*3);
#endif

    size_t _i=0;
    double dx=0,dy=0,dz=0;
    //iterate over all grids
    if(verbose)
      std::cout<<"Sampling deformations... "<<std::flush;
    
    
    for(size_t s=0;s<_mask.size();s++)
    {
      if(verbose)
        std::cout<<std::endl<<s<<" of "<<_mask.size()<<" "<<std::flush;
      minc::mask3d_iterator it(_mask[s], _mask[s]->GetRequestedRegion() );

      for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
        tag_point p1,p2,p3;

        if(!it.Value()) continue;

        minc::mask3d::IndexType idx=it.GetIndex();
        if(  skip_voxels>1 && 
            (idx[0]%skip_voxels || idx[1]%skip_voxels || idx[2]%skip_voxels )) 
          continue;

        //TODO: find the original point here from previous iteration
        _mask[s]->TransformIndexToPhysicalPoint(idx,p1);

        p2=_lin_xfm[s]->TransformPoint(p1);
        
        fun_x.generate_basis(bas_x1,order,p1);
        fun_x.generate_basis(bas_x2,order,p2);

        //TODO: Add rotation matrix here, from PCA?
        
        p3=_nl_xfm[s]->TransformPoint(p2);
        
        for(size_t _k=0; _k<3;_k++) //coordinates
        {
#ifndef SEPARATE    
          bas_diff.assign((order-3)*3,0.0);
#endif
          for(size_t _j=0; _j<(order-3); _j++) {
#ifdef SEPARATE
            bas_diff[_j]=bas_x2[_j+3]-bas_x1[_j+3];
#else
            bas_diff[_j+(order-3)*_k]= bas_x2[_j+3]-bas_x1[_j+3];
#endif
          }
          //TODO: apply rotation matrix here
#ifdef SEPARATE
          pol[_k].accumulate(bas_diff, p3[_k]-p2[_k]);
#else
          pol.accumulate(bas_diff, p3[_k]-p2[_k]);
#endif
        }
        _i++;

        dx+=(p3[0]-p2[0])*(p3[0]-p2[0]);
        dy+=(p3[1]-p2[1])*(p3[1]-p2[1]);
        dz+=(p3[2]-p2[2])*(p3[2]-p2[2]);

        if(!(_i%1000) && verbose)
          std::cout<<"."<<std::flush;
      }
    }
    dx/=_i;
    dy/=_i;
    dz/=_i;
    
    if(verbose)
      std::cout<<"dx="<<sqrt(dx)<<" dy="<<sqrt(dy)<<" dz="<<sqrt(dz)<<" Done!"<<std::endl;
    
    //TODO 
    if(verbose)
      std::cout<<"Fitting:"<<std::flush;
    
#ifdef SEPARATE
    fitting_coeff _coeff(order-3);
    
    for(size_t _k=0; _k<3;_k++) //coordinates
    {
      _coeff[_k].resize(order-3);
      pol[_k].solve_unstable(_coeff,0.01,verbose);
      for(size_t _j=3; _j<order; _j++) {
        coeff[_k][_j]=_coeff[_j-3];
      }
    }
#else
    fitting_coeff _coeff((order-3)*3);
    pol.solve_unstable(_coeff);
    
    for(size_t _j=3; _j<order; _j++) {
      for(size_t _k=0;_k<3;_k++)
        coeff[_k][_j]=_coeff[_j-3+(order-3)*_k];
    }
#endif    
    //TODO: convert back to spherical harmonics coeff

    //linear part is identity
    coeff[0][1]=coeff[1][2]=coeff[2][0]=1.0;
    coeff[0][0]=coeff[0][2]=0;
    coeff[1][0]=coeff[1][1]=0;
    coeff[2][1]=coeff[2][2]=0;
  }
  
  void save_residuals(const char *prefix) {
    fitting_coeff _coeff((order-3)*3);
    
    for(size_t _j=3; _j<order; _j++) {
      for(size_t _k=0;_k<3;_k++)
        _coeff[_j-3+(order-3)*_k]=coeff[_k][_j];
    }
    
    basis_vector bas_x1(order);
    basis_vector bas_x2(order);
    
    basis_vector bas_diff((order-3)*3);

    minc::MNK_Gauss_Polinomial pol((order-3)*3);

    size_t _i=0;
    //iterate over all grids
    if(verbose)
      std::cout<<"Sampling deformations... "<<std::flush;
    
    for(size_t s=0;s<_mask.size();s++)
    {
      minc::image3d::Pointer residual_x=minc::image3d::New();
      minc::image3d::Pointer residual_y=minc::image3d::New();
      minc::image3d::Pointer residual_z=minc::image3d::New();
      
      allocate_same<minc::image3d,minc::mask3d>(residual_x,_mask[s]);
      allocate_same<minc::image3d,minc::mask3d>(residual_y,_mask[s]);
      allocate_same<minc::image3d,minc::mask3d>(residual_z,_mask[s]);
      
      minc::mask3d_iterator    it(_mask[s],     _mask[s]->GetRequestedRegion() );
      minc::image3d_iterator it_x(residual_x, residual_x->GetRequestedRegion() );
      minc::image3d_iterator it_y(residual_y, residual_y->GetRequestedRegion() );
      minc::image3d_iterator it_z(residual_z, residual_z->GetRequestedRegion() );

      for(it.GoToBegin(),it_x.GoToBegin(),it_y.GoToBegin(),it_z.GoToBegin(); !it.IsAtEnd(); 
          ++it,++it_x,++it_y,++it_z) {
        
        tag_point p1,p2,p3;
        it_x.Set(0.0);it_y.Set(0.0);it_z.Set(0.0);

        if(!it.Value()) continue;

        minc::mask3d::IndexType idx=it.GetIndex();
        if(skip_voxels>1 && ( idx[0]%skip_voxels || idx[1]%skip_voxels || idx[2]%skip_voxels ) ) 
          continue;

        //TODO: find the original point here from previous iteration
        _mask[s]->TransformIndexToPhysicalPoint(idx,p1);
        p2=_lin_xfm[s]->TransformPoint(p1);

        fun_x.generate_basis(bas_x1,order,p1);
        fun_x.generate_basis(bas_x2,order,p2);
        bas_diff.assign((order-3)*3,0.0);

        p3=_nl_xfm[s]->TransformPoint(p2);

        double res[3];
        for(size_t _k=0; _k<3;_k++) //coordinates
        {
          for(size_t _j=0; _j<(order-3); _j++) {
            bas_diff[_j+(order-3)*_k] = bas_x2[_j+3] - bas_x1[_j+3];
          }
          //TODO: apply rotation matrix here
          res[_k]=pol.fit(bas_diff,_coeff)-(p3[_k]-p2[_k]);
        }
        it_x.Set(res[0]);
        it_y.Set(res[1]);
        it_z.Set(res[2]);
        
        if(!(_i%1000) && verbose)
          std::cout<<"."<<std::flush;
      }
      char tmp[1024];
#ifdef HAVE_MINC4ITK
      minc::set_minc_storage_type(residual_x,NC_SHORT,true);
      minc::set_minc_storage_type(residual_y,NC_SHORT,true);
      minc::set_minc_storage_type(residual_z,NC_SHORT,true);
#else
      minc::set_minc_storage_type(residual_x,typeid(unsigned short).name());
      minc::set_minc_storage_type(residual_y,typeid(unsigned short).name());
      minc::set_minc_storage_type(residual_z,typeid(unsigned short).name());
#endif      
      sprintf(tmp,"%s_%lu_x.mnc",prefix,s);
      if(verbose)
        std::cout<<"Saving:"<<tmp<<std::endl;
      save_minc<minc::image3d>(tmp, residual_x);
      
      sprintf(tmp,"%s_%lu_y.mnc",prefix,s);
      if(verbose)
        std::cout<<"Saving:"<<tmp<<std::endl;
      save_minc<minc::image3d>(tmp, residual_y);

      sprintf(tmp,"%s_%lu_z.mnc",prefix,s);
      if(verbose)
        std::cout<<"Saving:"<<tmp<<std::endl;
      save_minc<minc::image3d>(tmp, residual_z);
      
    }
  }

  bool fit_tags(bool condition=false) {
    _last_distance=0.0;
    int i=0;
    sd=0;
//      do {
    fit_coeff(condition);
//        i++;
    //std::cout<<"\t"<<i;
//      } while(i<_max_iterations && keep<1.0&&remove_outliers());
    return true;
  }

  void load_sample(const char *lin_xfm_f,const char *nl_xfm_f,const char *mask_f) {
    if(verbose)
      std::cout<<"Loading lin xfm:"<<lin_xfm_f<<", nl_xfm:"<<nl_xfm_f<<", mask:"<<mask_f<<
        " skipping:"<<skip_voxels<<" ... "<<std::flush;
    
    minc::mask3d::Pointer mask(minc::mask3d::New());
    XFMTransformType::Pointer lin_xfm(XFMTransformType::New());
    XFMTransformType::Pointer nl_xfm(XFMTransformType::New());
    
    load_minc<minc::mask3d>(mask_f, mask);
    lin_xfm->OpenXfm(lin_xfm_f);
    nl_xfm->OpenXfm(nl_xfm_f);
    
    _lin_xfm.push_back(lin_xfm);
    _nl_xfm.push_back(nl_xfm);
    _mask.push_back(mask);
    
    minc::mask3d_iterator it(mask, mask->GetRequestedRegion() );
    
    size_t cnt=0;

    for(it.GoToBegin(); !it.IsAtEnd(); ++it) {
      tag_point p;
      minc::mask3d::IndexType idx=it.GetIndex();
      mask->TransformIndexToPhysicalPoint(idx,p);
      minc::mask3d::IndexType idx_m;

      if(!it.Value())
        continue;

      if( skip_voxels>1 && 
        ( idx[0]%skip_voxels || idx[1]%skip_voxels || idx[2]%skip_voxels )) continue;

      cnt++;
    }
    _sample_count+=cnt;

    if(verbose)
      std::cout<<cnt<<" nodes"<<std::endl;
  }
};

const double Tag_fit::_distance_epsilon=1e-10;

void show_usage (const char * prog)
{
  std::cerr
      << "Usage: "<<prog<<" <lin_xfm 1> <nl_xfm 1> <mask 1> [<lin_xfm 2> <nl_xfm 2> <mask 2>.... <lin_xfm n> <nl_xfm n> <mask n>]" << endl
      << "--clobber overwrite files"    << endl
      << "--order <n> (3)"<<endl
      << "--keep <part> 0.0-1.0 part of data points to keep (0.8)"<<endl
      << "--iter <n> maximum number of iterations (200)"<<endl
      << "--skip <n> subsample deformation field by factor n"<<endl
      << "--lambda <f> set lambda for lasso solver, default 0.5"<<endl
      << "--output <output.par> , default output.par "<<endl
      << "--residual <base> , output residuals "<<endl
      << "--scaling <d> set scaling for Spherial Harmonics"<<endl;
}

int main (int argc, char **argv)
{
  int verbose=0, clobber=0,skip_grid=0,cond=0;
  double max=5.0;
  double keep=1.0;
  double remove=0;
  int iter=1;
  int order=3;
  int skip_voxels=0;
  std::string grid_f,mask_f,output_f,dump_f;
  std::string residual_f;
  int lsq=12;
  int limit_linear=0;
  double lambda=0.5;
  double scaling=200.0;
  

  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"order",   required_argument,   0, 'o'},
    {"keep",    required_argument,   0, 'k'},
    {"iter",    required_argument,   0, 'i'},
    {"skip",    required_argument,   0, 's'},
    {"version", no_argument,         0, 'v'},
    {"lambda",  required_argument,   0, 'l'},
    {"output",  required_argument,   0, 'O'},
    {"residual",  required_argument, 0, 'r'},
    {"scaling",    required_argument,   0, 'S'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "o:k:i:vs:l:O:r:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c) {
    case 0:
      break;
    case 'v':
      cout << "Version: 0.2" << endl;
      return 0;
    case 'o':
      order=atoi(optarg);
      break;
    case 'k':
      keep=atof(optarg);
      break;
    case 's':
      skip_voxels=atoi(optarg);
      break;
    case 'i':
      iter=atoi(optarg);
      break;
    case 'l':
      lambda=atof(optarg);
      break;
    case 'O':
      output_f=optarg;
      break;
    case 'r':
      residual_f=optarg;
      break;
    case 'S':
      scaling=atof(optarg);break;
    case '?':
      /* getopt_long already printed an error message. */
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if ((argc - optind) < 3 || (argc - optind)%3 ) {
    show_usage (argv[0]);
    return 1;
  }
  if(output_f.empty())
    output_f="output.par";

  try {
    float max_dev=10;
    Tag_fit fit(basis_functions_x::parameters_no(order),keep,max_dev,iter,verbose,false,skip_voxels);

    fit.scaling=scaling;
    fit.lambda=lambda;

    if (!clobber && !access(output_f.c_str (), F_OK)) {
      cerr << output_f.c_str () << " Exists!" << endl;
      return 1;
    }
    
    while((argc - optind)>0) {
      fit.load_sample(argv[optind],argv[optind+1],argv[optind+2]);
      optind+=3;
    }

    if(!fit.fit_tags(cond)) {
      cerr<<"Fitting failed, due to large Standard Deviation !"<<endl;
      return 10;
    }

    if(!residual_f.empty())
      fit.save_residuals(residual_f.c_str());
    
    std::ofstream cf(output_f.c_str());

    if(!cf.good())
      REPORT_ERROR("Can't open file for writing!");

    save_coeff(cf,fit.coeff[0]);
    save_coeff(cf,fit.coeff[1]);
    save_coeff(cf,fit.coeff[2]);

  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  } catch( itk::ExceptionObject & err ) {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
  return 0;

}
