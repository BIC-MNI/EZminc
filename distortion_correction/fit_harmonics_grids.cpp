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

#ifdef HAVE_MINC4ITK
#include <itkMincHelpers.h>
#else
#include "itk4MincHelpers.h"
#endif

using namespace std;
using namespace minc;

void print_coeff (std::ostream & out, const std::vector < double >&coeff2)
{
	out.precision(40);
	int order = coeff2.size ();
	for (int i = 0; i < order; i++)
	{
		out << coeff2[i] << " ";
		if(!out.good())
			ITK_REPORT_ERROR("Can't write to file");
		
	}
	out << std::endl;
}

void save_coeff (std::ostream & out, const std::vector < double >&coeff2)
{
	int order = coeff2.size ();
	out.precision(40);
	for (int i = 0; i < order; i++)
	{
		out << coeff2[i] << endl;
		if(!out.good())
			ITK_REPORT_ERROR("Can't write to file");
		
	}
}

class Tag_fit
{
	protected:
    double _last_distance;
    static const double _distance_epsilon;
  public:
    typedef std::vector <bool> fitting_mask;
    typedef std::vector <double> fitting_coeff;
    typedef std::vector <fitting_coeff> fittings;
    typedef std::vector <int> Index;
    typedef std::vector <double> Distances;
    
    int       order;
    double    keep;
    bool      verbose;
    double    max_dev;
    int       _max_iterations;
    double    sd,max_distance;
    fittings  coeff;
    bool      limit_linear;
    bool      cache_basis;
    int       skip_voxels;
    tag_points ideal, measured;
    fitting_mask mask;
    fittings  basis_x,basis_y,basis_z;
    Index     index;
    Distances distances;

		basis_functions_x fun_x; // here we have the same basis for X,Y,Z


    Tag_fit(int order_, double keep_, double max_dev_=5.0, int max_iter=100,bool verbose_=false,bool ll=false,int skip_voxels_=0):
           order(order_), 
           keep(keep_), 
           verbose(verbose_), 
           max_dev(max_dev_),
           _max_iterations(max_iter),
           sd(0), 
           max_distance(0), 
           coeff(3),
           limit_linear(ll),
           cache_basis(max_iter>1),
           skip_voxels(skip_voxels_)
    {
    }


    void reset_index(void)
    {
      index.resize(ideal.size());
      for(size_t i=0;i<ideal.size();i++)
        index[i]=i;
    }
    
    void calculate_basis(void)
    {
      if(!cache_basis) return;
      basis_x.resize(ideal.size());
      //basis_y.resize(ideal.size());
      //basis_z.resize(ideal.size());
      tag_points::const_iterator i=ideal.begin();
      
      fittings::iterator mx=basis_x.begin();
      //fittings::iterator my=basis_y.begin();
      //fittings::iterator mz=basis_z.begin();
      basis_vector bas(order);
      for(; i!=ideal.end(); i++, mx++/*,my++,mz++*/)
      {
        mx->resize(order);
        //my->resize(order);
        //mz->resize(order);
        
        //all basis are the same
        fun_x.generate_basis(bas,order,*i);
        for(int k=0;k<order;k++)
        {
          (*mx)[k]=bas[k];
        }
      }
      
    }
    
    
    
    void fit_coeff(bool condition=false)
    {
      if(ideal.size()!=measured.size())
        ITK_REPORT_ERROR("Mismatching number of points!");
      if(ideal.size()!=mask.size())
        ITK_REPORT_ERROR("Mismatching number of points!");
    
      reset_index();
      coeff.resize(3);
      coeff[0].resize(order);
      coeff[1].resize(order);
      coeff[2].resize(order);
      
      minc::MNK_Gauss_Polinomial pol_x(limit_linear?order-3:order);
      minc::MNK_Gauss_Polinomial pol_y(limit_linear?order-3:order);
      minc::MNK_Gauss_Polinomial pol_z(limit_linear?order-3:order);
      
      tag_points::const_iterator j=measured.begin();
      tag_points::const_iterator i=ideal.begin();
      fitting_mask::const_iterator m=mask.begin();
      
      fittings::const_iterator bx=basis_x.begin();
      basis_vector bas_x(order);
      
      for(; i!=ideal.end(); i++, j++, m++)
      {
        if(!*m)
        {
          if(cache_basis)
          {
            bas_x=*bx;
          } else {
            fun_x.generate_basis(bas_x,order,*i);
          }
          
          if(limit_linear)
          { 
            //TODO: fix indexing
            basis_vector _bas(bas_x.begin()+3,bas_x.end());
            
            pol_x.accumulate(_bas, (*j)[0]-(*i)[0]);
            pol_y.accumulate(_bas, (*j)[1]-(*i)[1]);
            pol_z.accumulate(_bas, (*j)[2]-(*i)[2]);
          } else {
            pol_x.accumulate(bas_x, (*j)[0]);
            pol_y.accumulate(bas_x, (*j)[1]);
            pol_z.accumulate(bas_x, (*j)[2]);
          }
        }

        if(cache_basis)
        {
          bx++; 
        }
      }
      
      double cond_x,cond_y,cond_z;
      if(limit_linear)
      {
        fittings  _coeff(3);
        
        _coeff[0].resize(order-3);
        _coeff[1].resize(order-3);
        _coeff[2].resize(order-3);
        
        cond_x=pol_x.solve_unstable(_coeff[0],0.01,verbose);
        cond_y=pol_y.solve_unstable(_coeff[1],0.01,verbose);
        cond_z=pol_z.solve_unstable(_coeff[2],0.01,verbose);
        
        coeff[0][1]=coeff[1][2]=coeff[2][0]=1.0;
        coeff[0][0]=coeff[0][2]=0;
        coeff[1][0]=coeff[1][1]=0;
        coeff[2][1]=coeff[2][2]=0;
        
        for(size_t _j=3; _j<order; _j++) {
          for(size_t _k=0;_k<3;_k++)
            coeff[_k][_j]=_coeff[_k][_j-3];
        }
        
      } else {
        cond_x=pol_x.solve_unstable(coeff[0],0.01,verbose);
        cond_y=pol_y.solve_unstable(coeff[1],0.01,verbose);
        cond_z=pol_z.solve_unstable(coeff[2],0.01,verbose);
      }
      
      if(condition)
      {
        cout<<"cond_x="<<cond_x<<"\t";
        cout<<"cond_y="<<cond_y<<"\t";
        cout<<"cond_z="<<cond_z<<endl;
      }
    }
	
    class IndexSort
    {
      Distances &distances;
      public:
      IndexSort(Distances &distances_):
      distances(distances_)
      {}
      bool operator()(int i,int j)
      {
        return distances[i]<distances[j];
      }
    };
    
    void build_index(void)
    {
      IndexSort _sort(distances);
      std::sort<Index::iterator,IndexSort>(index.begin(), index.end(), _sort);
    }
    
    double evaluate_distance(double keep)
    {
      int size=ideal.size()*keep;
      double sd=0.0;
      for(int i=0;i<ideal.size();i++)
      {
        if(i<size) {
          sd+=distances[index[i]];
          mask[index[i]]=false;
        } else mask[index[i]]=true;
      }
      sd/=size;
      return sqrt(sd);
    }
    
    bool remove_outliers(void)
    {
      minc::MNK_Gauss_Polinomial pol_x(order);
      minc::MNK_Gauss_Polinomial pol_y(order);
      minc::MNK_Gauss_Polinomial pol_z(order);
      
      distances.resize(ideal.size());
      tag_points::const_iterator j=measured.begin();
      tag_points::const_iterator i=ideal.begin();
      
      fittings::const_iterator bx=basis_x.begin();
      int k=0;
      int max_k=-1;
      int cnt=0;
      max_distance=0.0;
      sd=0.0;
      basis_vector bas_x(order);

      for(;i!=ideal.end();i++, j++, k++)
      {
        if(cache_basis)
        {
          bas_x=*bx;
        } else {
          fun_x.generate_basis(bas_x,order,*i);
        }
        
        cnt++;
        tag_point moved;
        moved[0]=pol_x.fit(bas_x, coeff[0]);
        moved[1]=pol_y.fit(bas_x, coeff[1]);
        moved[2]=pol_z.fit(bas_x, coeff[2]);
        
        double d=(*j).SquaredEuclideanDistanceTo(moved);
        distances[k]=d;
        
        if(!mask[k]&&d>max_distance) { max_distance=d; max_k=k;}
        
        if(cache_basis)
        {
          bx++;
        }        
      }
      max_distance=sqrt(max_distance);
      build_index();
      sd=evaluate_distance(keep);
      if(verbose)
        cout<<sd<<":"<<max_distance<<"\t";

      return true;
    }
    
    bool fit_tags(bool condition=false)
    {
      _last_distance=0.0;
      if(measured.size()!=ideal.size()) 
        ITK_REPORT_ERROR("Mismatching number of tags!");
      mask.resize(ideal.size(),false);
      calculate_basis();
      int i=0;
      sd=0;
      do {
        fit_coeff(condition);
        i++;
        //std::cout<<"\t"<<i;
      } while(i<_max_iterations && keep<1.0&&remove_outliers());
      return sd<max_dev;
    }
    
    void load_grid(const char *grid_f,const char *mask_f)
    {
      if(verbose)
        std::cout<<"Loading grid:"<<grid_f<<" and mask:"<<mask_f<<" ... ";
      
      minc::def3d::Pointer grid (minc::def3d::New());
      minc::mask3d::Pointer mask(minc::mask3d::New());
      
      load_minc<minc::def3d>(grid_f, grid);
      load_minc<minc::mask3d>(mask_f, mask);
      
      minc::def3d_iterator it(grid, grid->GetRequestedRegion() );
      int cnt=0;

      for(it.GoToBegin();!it.IsAtEnd();++it)
      {
        tag_point p;
        minc::def3d::IndexType idx=it.GetIndex();
        grid->TransformIndexToPhysicalPoint(idx,p);
        minc::mask3d::IndexType idx_m;
        
        if(skip_voxels>1 && ( idx[0]%skip_voxels || idx[1]%skip_voxels || idx[2]%skip_voxels )) continue;
        
        if(!mask->TransformPhysicalPointToIndex(p,idx_m) || !mask->GetPixel(idx_m)) 
          continue;
        
        
        
        ideal.push_back(p);
        p[0]+=it.Value()[0];
        p[1]+=it.Value()[1];
        p[2]+=it.Value()[2];
        measured.push_back(p);
        cnt++;
      }
      if(verbose)
        std::cout<<cnt<<" nodes"<<std::endl;
    }
};

const double Tag_fit::_distance_epsilon=1e-10;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" <grid1> <mask1> [<grid2> <mask2> .... <grid n> <mask n>] <output.par> " << endl
    << "--clobber overwrite files"    << endl
    << "--order <n> (3)"<<endl
    << "--keep <part> 0.0-1.0 part of data points to keep (0.8)"<<endl
    << "--iter <n> maximum number of iterations (200)"<<endl
    << "--remove <pct> 0-1 (0) randomly remove voxels"<<endl
    << "--cond calculate condition number"<<endl
    << "--limit limit linear component to identity"<<endl
    << "--skip <n> subsample deformation field by factor n"<<endl;
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
  std::string grid_f,mask_f,output,dump_f;
  std::string residuals_f;
  int lsq=12;
  int limit_linear=0;
  static struct option long_options[] = {
		{"verbose", no_argument,       &verbose, 1},
		{"quiet",   no_argument,       &verbose, 0},
		{"clobber", no_argument,       &clobber, 1},
    {"order",   required_argument,   0, 'o'},
    {"keep",    required_argument,   0, 'k'},
    {"iter",    required_argument,   0, 'i'},
    {"skip",    required_argument,   0, 's'},
		{"version", no_argument,         0, 'v'},
    {"cond", no_argument,       &cond, 1},
    {"limit", no_argument,       &limit_linear, 1},
    //{"remove",  required_argument,  0, 'e'},
		{0, 0, 0, 0}
		};
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "o:k:i:vs:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
			{
			case 0:
				break;
			case 'v':
				cout << "Version: 0.2" << endl;
				return 0;
      case 'o':
        order=atoi(optarg);break;
      case 'k':
        keep=atof(optarg);break;
      case 'i':
        iter=atoi(optarg);break;
      case 's':
        skip_voxels=atoi(optarg);break;
			case '?':
				/* getopt_long already printed an error message. */
			default:
				show_usage (argv[0]);
				return 1;
			}
    }

	if ((argc - optind) < 3 || !((argc - optind)&1) ) {
		show_usage (argv[0]);
		return 1;
	}
	try
  {
    gsl_rng_env_setup();
    
    float max_dev=10;
    Tag_fit fit(basis_functions_x::parameters_no(order),keep,max_dev,iter,verbose,limit_linear,skip_voxels);
  
    while((argc - optind)>1)
    {
      fit.load_grid(argv[optind],argv[optind+1]);
      optind+=2;
    }
    output=argv[optind];
    if (!clobber && !access(output.c_str (), F_OK))
    {
      cerr << output.c_str () << " Exists!" << endl;
      return 1;
    }
    if(!fit.fit_tags(cond))
    {
      cerr<<"Fitting failed, due to large Standard Deviation !"<<endl;
      return 10;
    }
    
    std::ofstream cf(output.c_str());
    
    if(!cf.good())
      ITK_REPORT_ERROR("Can't open file for writing!");
    
    save_coeff(cf,fit.coeff[0]);
    save_coeff(cf,fit.coeff[1]);
    save_coeff(cf,fit.coeff[2]);
    
  }
#ifdef HAVE_MINC4ITK
  catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1; 
  }
#endif  
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
	return 0;
	
}
