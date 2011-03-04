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
#include "minc_wrappers.h"
//#include "minc_io.h"
#include <iostream>
#include <fstream>
#include "data_proc.h"
#include <itkBSplineInterpolateImageFunction.h>
#include <unistd.h>
#include <getopt.h>
#include "sphericalHarmonicsTransform.h"
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>

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
			REPORT_ERROR("Can't write to file");
		
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
			REPORT_ERROR("Can't write to file");
		
	}
	//out << std::endl;
}

class Tag_fit
{
	protected:
    double _last_distance;
    static const double _distance_epsilon;
    public:
    int _max_iterations;
    typedef std::vector <bool> fitting_mask;
    typedef std::vector <double> fitting_coeff;
    typedef std::vector <fitting_coeff> fittings;
    typedef std::vector <int> Index;
    typedef std::vector <double> Distances;
    
    
    tag_points ideal, measured;
    fitting_mask mask;
    bool     verbose;
    bool     limit_linear;
    double   keep;
    double   max_dev;
    double   sd,max_distance;
    int      order;
    fittings coeff;
    fittings basis_x;//,basis_z;
    Index    index;
    Distances distances;

		CylindricalFunctions fun_x;
	
    void reset_index(void)
    {
      index.resize(ideal.size());
      for(int i=0;i<ideal.size();i++)
        index[i]=i;
    }
    
    void calculate_basis(void)
    {
      basis_x.resize(ideal.size());
      tag_points::const_iterator i=ideal.begin();
      
      fittings::iterator mx=basis_x.begin();
      basis_vector bas(order);
      for(; i!=ideal.end(); i++, mx++)
      {
        mx->resize(order);
        
        //all basis are the same
        fun_x.generate_basis(bas,order,*i);
        for(int k=0;k<order;k++)
        {
          (*mx)[k]=bas[k];
        }
      }
      
    }
    
    Tag_fit(int order_, double keep_, double max_dev_=5.0, int max_iter=100,bool verbose_=false,bool ll=false):
      order(order_), keep(keep_), verbose(verbose_), max_dev(max_dev_),_max_iterations(max_iter),
      max_distance(0), sd(0), coeff(3),limit_linear(ll)
    {
    }
    
    
    void fit_coeff(bool condition=false)
    {
      if(ideal.size()!=measured.size())
        REPORT_ERROR("Mismatching number of points!");
      if(ideal.size()!=mask.size())
        REPORT_ERROR("Mismatching number of points!");
    
      reset_index();
      coeff.resize(2);
      coeff[0].resize(order);
      coeff[1].resize(order);
      
      minc::MNK_Gauss_opt<tag_point> pol_x(limit_linear?order-2:order);
      minc::MNK_Gauss_opt<tag_point> pol_z(limit_linear?order-2:order);
      
      tag_points::const_iterator j=measured.begin();
      tag_points::const_iterator i=ideal.begin();
      fitting_mask::const_iterator m=mask.begin();
      
      fittings::const_iterator bx=basis_x.begin();
      for(; i!=ideal.end(); i++, j++, m++, bx++)
      {
        if(*m) continue;
          //this is a quick hack
        if(limit_linear)
        {
          double dx=(*j)[0]-(*i)[0];
          double dy =(*j)[1]-(*i)[1];
          pol_x.accumulate(&(*bx)[2], sqrt(dx*dx+dy*dy));
          pol_z.accumulate(&(*bx)[2], (*j)[2]-(*i)[2]);
        } else {
          pol_x.accumulate(*bx, sqrt((*j)[0]*(*j)[0]+(*j)[1]*(*j)[1]));
          pol_z.accumulate(*bx, (*j)[2]);
        }
      }
      //if(condition)
      //{
      double cond_x,cond_y,cond_z;
      if(limit_linear)
      {
        cond_x=pol_x.solve_svd(&coeff[0][2]);
        cond_z=pol_z.solve_svd(&coeff[1][2]);
        
        coeff[0][0]=coeff[1][1]=1.0;
        coeff[0][1]=coeff[1][0]=0;
      } else {
        cond_x=pol_x.solve_svd(coeff[0]);
        cond_z=pol_z.solve_svd(coeff[1]);
      }
      
      if(condition)
      {
        cout<<"cond_x="<<cond_x<<"\t";
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
      minc::MNK_Gauss < tag_point, basis_functions_x> pol_x(order);
      minc::MNK_Gauss < tag_point, basis_functions_z> pol_z(order);
      distances.resize(ideal.size());
      tag_points::const_iterator j=measured.begin();
      tag_points::const_iterator i=ideal.begin();
      fittings::const_iterator bx=basis_x.begin();
      //fittings::const_iterator bz=basis_z.begin();
      int k=0;
      int max_k=-1;
      int cnt=0;
      max_distance=0.0;
      double dx=0.0,dy=0.0,dz=0.0;
      sd=0.0;
      for(;i!=ideal.end();i++, j++, k++, bx++/*, bz++*/)
      {
        //if(mask[k]) continue;
        cnt++;
        tag_point moved;
        double r=pol_x.fit(*bx, coeff[0], *i);
        double ir=sqrt((*i)[0]*(*i)[0]+(*i)[1]*(*i)[1]);
        if(ir>1e-6)
        {
          moved[0]=(*i)[0]*r/ir;
          moved[1]=(*i)[1]*r/ir;
        } else {
          moved[0]=(*i)[0];
          moved[1]=(*i)[1];
        }
        moved[2]=pol_z.fit(*bx, coeff[1], *i);
        double sx=(*j)[0]-moved[0];
        double sy=(*j)[1]-moved[1];
        double sz=(*j)[2]-moved[2];
        dx+=sx*sx;
        dy+=sy*sy;
        dz+=sz*sz;
        double d=sx*sx+sy*sy+sz*sz;//
        distances[k]=d;
        //sd+=d;
        if(!mask[k]&&d>max_distance) { max_distance=d; max_k=k;}
      }
      if(cnt>0)
      {
        dx=sqrt(dx/cnt);
        dy=sqrt(dy/cnt);
        dz=sqrt(dz/cnt);
      }
      
      max_distance=sqrt(max_distance);
      build_index();
      sd=evaluate_distance(keep);
      if(verbose)
        cout<<sd<<","<<max_distance<<","<<dx<<","<<dy<<","<<dz<<std::endl;
      return true;
    }
    
    bool fit_tags(bool condition=false)
    {
      _last_distance=0.0;
      if(measured.size()!=ideal.size()) 
        REPORT_ERROR("Mismatching number of tags!");
      mask.resize(ideal.size(),false);
      calculate_basis();
      int i=0;
      do {
        fit_coeff(condition);
        i++;
        //std::cout<<"\t"<<i;
      } while(remove_outliers() && i<_max_iterations && keep<1.0);
      return sd<max_dev;
    }
    
    void load_grid(const char *grid_f,const char *mask_f)
    {
      if(verbose)
        std::cout<<"Loading grid:"<<grid_f<<" and mask:"<<mask_f<<" ... ";
      minc::def3d::Pointer grid (minc::def3d::New());
      minc::mask3d::Pointer mask(minc::mask3d::New());
      
      load_minc(grid_f, grid);
      load_minc(mask_f, mask);
      
      minc::def3d_iterator it(grid, grid->GetRequestedRegion() );
      int cnt=0;

      for(it.GoToBegin();!it.IsAtEnd();++it)
      {
        tag_point p;
        minc::def3d::IndexType idx=it.GetIndex();
        grid->TransformIndexToPhysicalPoint(idx,p);
        minc::mask3d::IndexType idx_m;
        if(!mask->TransformPhysicalPointToIndex(p,idx_m) || !mask->GetPixel(idx_m)) continue;

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
    << "--limit limit linear component to identity"<<endl;
    //<< "--scale <d>"<<endl;
}

int main (int argc, char **argv)
{
  int verbose=0, clobber=0,skip_grid=0,cond=0;
  double max=5.0;
  double keep=0.8;
  double remove=0;
  double scale=100;
  int iter=200;
  int order=3;
  std::string grid_f,mask_f,output,dump_f;
  std::string residuals_f;
  //TODO: fix this to a proper default
  int lsq=12;
  int limit_linear=0;
  static struct option long_options[] = {
		{"verbose", no_argument,       &verbose, 1},
		{"quiet",   no_argument,       &verbose, 0},
		{"clobber", no_argument,       &clobber, 1},
    {"order",   required_argument,   0, 'o'},
    {"keep",    required_argument,   0, 'k'},
    {"iter",    required_argument,   0, 'i'},
		{"version", no_argument,         0, 'v'},
    {"cond", no_argument,       &cond, 1},
    {"limit", no_argument,       &limit_linear, 1},
    //{"scale", required_argument,      0, 's'},
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
      case 's':
        scale=atof(optarg);break;
      case 'i':
        iter=atoi(optarg);break;
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
		Tag_fit fit(CylindricalFunctions::parameters_no(order),keep,max_dev,iter,verbose,limit_linear);
  
		minc::def3d::Pointer grid(minc::def3d::New());
		minc::mask3d::Pointer mask(minc::mask3d::New());
    
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
      REPORT_ERROR("Can't open file for writing!");
    save_coeff(cf,fit.coeff[0]);
    save_coeff(cf,fit.coeff[1]);
    
	} catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1; 
  }
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
	return 0;
	
}
