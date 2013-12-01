/* ----------------------------- MNI Header -----------------------------------
@NAME       :  fit_harmonics_grids_regularize
@DESCRIPTION:  spherical harminic LSQ approximation programm with legendre regularization
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

#ifdef HAVE_MINC4ITK
#include <itkMincHelpers.h>
#else
#include "itk4MincHelpers.h"
#endif

#include "sphericalHarmonicsTransform.h"

#include <gsl/gsl_rng.h>
#include <time.h>
#include <sys/time.h>

using namespace std;
using namespace minc;

void print_coeff(std::ostream & out, const std::vector < double >&coeff)
{
	out.precision(40);
	int order = coeff.size ();
	for (int i = 0; i < order; i++)
	{
		out << coeff[i] << " ";
		if(!out.good())
			REPORT_ERROR("Can't write to file");
		
	}
	out << std::endl;
}

void save_coeff(std::ostream & out, const std::vector < double >&coeff)
{
	int order = coeff.size ();
	out.precision(40);
	for (int i = 0; i < order; i++)
	{
		out << coeff[i] << endl;
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
    typedef std::vector <bool> fitting_mask;
    typedef std::vector <double> fitting_coeff;
    typedef std::vector <fitting_coeff> fittings;
    typedef std::vector <int> Index;
    typedef std::vector <double> Distances;
    std::vector< std::vector<double> > pca_matrix;
    
    int      order,order_c;
    int _max_iterations;
    bool     verbose;
    bool     cylindric;
    double   max_dev;
    double   sd,max_distance;
    int      pcs;
    double legendre_coeff;
    int skip_voxels;
    fitting_coeff coeff;
    fittings basis_x;
    tag_points ideal, measured;
    fitting_mask mask;
    basis_vector regularize;
    basis_vector regularize_c;
    Index    index;
    Distances distances;

		basis_functions_x    fun_x;
    CylindricalFunctions fun_c;
	
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
      
      if(cylindric)
      {
        basis_vector bas(order_c);
        for(; i!=ideal.end(); i++, mx++)
        {
          mx->resize(order);
          
          fun_c.generate_basis(bas,order_c,*i);
          for(int k=0;k<order;k++)
          {
            (*mx)[k]=bas[k];
          }
        }
      } else {
        basis_vector bas(order);
        for(; i!=ideal.end(); i++, mx++)
        {
          mx->resize(order);
          
          fun_x.generate_basis(bas,order,*i);
          for(int k=0;k<order;k++)
          {
            (*mx)[k]=bas[k];
          }
        }
      }
      //used for regularization
      regularize.resize(order);
      regularize_c.resize(order_c);
      
      if(legendre_coeff>0.0)
      {
        fun_x.generate_regularization_vector(regularize,order,legendre_coeff);
        fun_c.generate_regularization_vector(regularize_c,order_c,legendre_coeff);
        
      }
      std::cout<<std::endl<<legendre_coeff<<std::endl;
    }
    
    Tag_fit(int order_, 
            int order_c_,
            double legendre,
            double max_dev_=10.0, 
            bool verbose_=false,
            bool cylindric_=false,
            int skip_voxels_=0
           ):
           order(order_),
           order_c(order_c_), 
           verbose(verbose_), 
           max_dev(max_dev_),
           max_distance(0), 
           sd(0), 
          legendre_coeff(legendre),
          pcs(-1),
          cylindric(cylindric_),
          skip_voxels(skip_voxels_)
    {
    }
    
    void fit_coeff(bool condition=false)
    {
      if(ideal.size()!=measured.size())
        REPORT_ERROR("Mismatching number of points!");
      if(ideal.size()!=mask.size())
        REPORT_ERROR("Mismatching number of points!");
      reset_index();
      
      if(cylindric)
      {
        if(verbose)
          std::cout<<"Performing cylindrical fit"<<std::endl;
        
        if(pcs>0)
        {
          //going to use principal components
          coeff.resize(order_c*2);
          
          std::vector<double> tmp[2];
          std::vector<double> tmp2(pcs);
          std::vector<double> coeff_pcs(pcs);
          
          tmp[0].resize(order_c*2);
          tmp[1].resize(order_c*2);
          
          
          minc::MNK_Gauss_Polinomial pol_pcs(pcs);
          
          tag_points::const_iterator j=measured.begin();
          tag_points::const_iterator i=ideal.begin();
          fitting_mask::const_iterator m=mask.begin();
          fittings::const_iterator bx=basis_x.begin();
          
          for(; i!=ideal.end(); i++, j++, m++, bx++)
          {
            if(*m) continue;
              //this is a quick hack
            
            for(int k=0;k<2;k++)
            {
              tmp[k].assign(order_c*2,0.0);
              for(int t=0;t<order_c;t++)
                tmp[k][t+k*order_c]=(*bx)[t];
            }
            
            //calculcate PCS basis
            double cyl[2];
            cyl[0]=sqrt((*j)[0]*(*j)[0]+(*j)[1]*(*j)[1]);
            cyl[1]=(*j)[2];
            
            for(int k=0;k<2;k++)
            {
              for(int t=0;t<pcs;t++)
              {
                tmp2[t]=0;
                for(int e=0;e<(order_c*2);e++)
                  tmp2[t]+=tmp[k][e]*pca_matrix[e][t];
              }
              pol_pcs.accumulate(tmp2,cyl[k]);
            }
          }
          
          pol_pcs.solve_unstable(coeff_pcs,0.01,verbose);
          //convert back to original coeffecients
          std::cout<<"Solution:";
          for(int j=0;j<pcs;j++)
          {
            std::cout<<coeff_pcs[j]<<" ";
          }
          std::cout<<std::endl;
          
          for(int e=0;e<(order_c*2);e++)
          {
            coeff[e]=0;
            for(int j=0;j<pcs;j++)
            {
              coeff[e]+=coeff_pcs[j]*pca_matrix[e][j];
            }
          }
        } else {
          coeff.resize(order_c*2);
          std::vector<double> tmp(order_c*2);
          
          minc::MNK_Gauss_Polinomial pol_x(order_c*2);
          
          tag_points::const_iterator j=measured.begin();
          tag_points::const_iterator i=ideal.begin();
          fitting_mask::const_iterator m=mask.begin();
          
          fittings::const_iterator bx=basis_x.begin();
          
          for(; i!=ideal.end(); i++, j++, m++, bx++)
          {
            if(*m) continue;
              //this is a quick hack
            
            double cyl[2];
            cyl[0]=sqrt((*j)[0]*(*j)[0]+(*j)[1]*(*j)[1]);
            cyl[1]=(*j)[2];
            
            for(int k=0;k<2;k++)
            {
              tmp.assign(tmp.size(),0.0);
              for(int t=0;t<order_c;t++)
                tmp[t+k*order_c]=(*bx)[t];
              
              pol_x.accumulate(tmp,cyl[k]);
            }
          }
          
          //now add regularization coeffecients
          if(legendre_coeff>0.0)
          {
            for(int j=0;j<order_c;j++)
            {
              pol_x.alpha().set(j,j,pol_x.alpha().get(j,j)+regularize_c[j]);
              pol_x.alpha().set(j+order_c,j+order_c,pol_x.alpha().get(j+order_c,j+order_c)+regularize_c[j]);
            }
          }
          pol_x.solve_unstable(coeff,0.01,verbose);
        }
      } else { //non-cylindrical case
      
        if(pcs>0)
        {
          //going to use principal components
          coeff.resize(order*3);
          
          std::vector<double> tmp[3];
          std::vector<double> tmp2(pcs);
          std::vector<double> coeff_pcs(pcs);
          
          tmp[0].resize(order*3);
          tmp[1].resize(order*3);
          tmp[2].resize(order*3);
          
          
          minc::MNK_Gauss_Polinomial pol_pcs(pcs);
          
          tag_points::const_iterator j=measured.begin();
          tag_points::const_iterator i=ideal.begin();
          fitting_mask::const_iterator m=mask.begin();
          fittings::const_iterator bx=basis_x.begin();
          
          for(; i!=ideal.end(); i++, j++, m++, bx++)
          {
            if(*m) continue;
              //this is a quick hack
            
            for(int k=0;k<3;k++)
            {
              tmp[k].assign(order*3,0.0);
              for(int t=0;t<order;t++)
                tmp[k][t+k*order]=(*bx)[t];
            }
            
            //calculcate PCS basis
            for(int k=0;k<3;k++)
            {
              for(int t=0;t<pcs;t++)
              {
                tmp2[t]=0;
                for(int e=0;e<(order*3);e++)
                  tmp2[t]+=tmp[k][e]*pca_matrix[e][t];
              }
              pol_pcs.accumulate(tmp2,(*j)[k]);
            }
          }
          pol_pcs.solve_unstable(coeff_pcs,0.01,verbose);
          //convert back to original coeffecients
          std::cout<<"Solution:";
          for(int j=0;j<pcs;j++)
          {
            std::cout<<coeff_pcs[j]<<" ";
          }
          std::cout<<std::endl;
          
          for(int e=0;e<(order*3);e++)
          {
            coeff[e]=0;
            for(int j=0;j<pcs;j++)
            {
                coeff[e]+=coeff_pcs[j]*pca_matrix[e][j];
            }
          }
        } else {
          coeff.resize(order*3);
          std::vector<double> tmp(order*3);
          
          minc::MNK_Gauss_Polinomial pol_x(order*3);
          
          tag_points::const_iterator j=measured.begin();
          tag_points::const_iterator i=ideal.begin();
          fitting_mask::const_iterator m=mask.begin();
          
          fittings::const_iterator bx=basis_x.begin();
          
          for(; i!=ideal.end(); i++, j++, m++, bx++/*, by++, bz++*/)
          {
            if(*m) continue;
            
            for(int k=0;k<3;k++)
            {
              tmp.assign(tmp.size(),0.0);
              for(int t=0;t<order;t++)
                tmp[t+k*order]=(*bx)[t];
              
              pol_x.accumulate(tmp,(*j)[k]);
            }
            
          }
          
          //now add regularization coeffecients
          if(legendre_coeff>0.0)
          {
            for(int j=0;j<order;j++)
            {
              pol_x.alpha().set(j,j,pol_x.alpha().get(j,j)+regularize[j]);
              pol_x.alpha().set(j+order,j+order,pol_x.alpha().get(j+order,j+order)+regularize[j]);
              pol_x.alpha().set(j+order*2,j+order*2,pol_x.alpha().get(j+order*2,j+order*2)+regularize[j]);
            }
          }
          pol_x.solve(coeff);
        }
      } 
      
    }
	
    bool fit_tags(bool condition=false)
    {
      _last_distance=0.0;
      if(measured.size()!=ideal.size()) 
        REPORT_ERROR("Mismatching number of tags!");
      mask.resize(ideal.size(),false);
      calculate_basis();
      int i=0;
      fit_coeff(condition);
      return sd<max_dev;
    }
    
    //TODO: works only for the first image!!!
    void save_error(const char *err_f,const char *grid_f)
    {
      
      minc::def3d::Pointer   grid=minc::load_minc<minc::def3d>(grid_f);
      minc::image3d::Pointer error(minc::image3d::New());
      allocate_same<minc::image3d,minc::def3d>(error,grid);
      error->FillBuffer(0.0);
      
      if(cylindric)
      {
        minc::MNK_Gauss_Polinomial pol_x(order_c*2);
        
        
        tag_points::const_iterator j=measured.begin();
        tag_points::const_iterator i=ideal.begin();
  
        fittings::const_iterator bx=basis_x.begin();
        
        std::vector<double> tmp(order_c*3);
        
        for(;i!=ideal.end();i++, j++, bx++)
        {
          
          minc::image3d::IndexType idx;
          error->TransformPhysicalPointToIndex(*i,idx);
          
          tag_point moved;
          double cyl[2];
          
          for(int k=0;k<2;k++)
          {
            tmp.assign(tmp.size(),0.0);
            
            for(int t=0;t<order_c;t++)
              tmp[t+k*order_c]=(*bx)[t];
            
            cyl[k]=pol_x.fit(tmp,coeff);
          }
          
          double ir=sqrt((*i)[0]*(*i)[0]+(*i)[1]*(*i)[1]);
          
          if(ir>1e-6)
          {
            moved[0]=(*i)[0]*cyl[0]/ir;
            moved[1]=(*i)[1]*cyl[0]/ir;
          } else {
            moved[0]=(*i)[0];
            moved[1]=(*i)[1];
          }
          moved[2]=cyl[1];
          
          error->SetPixel(idx,(*j).EuclideanDistanceTo(moved));
        }
      } else {
        minc::MNK_Gauss_Polinomial pol_x(order*3);
        
        minc::def3d::Pointer   grid=minc::load_minc<minc::def3d>(grid_f);
        minc::image3d::Pointer error(minc::image3d::New());
        
        
        error->FillBuffer(0.0);
        
        tag_points::const_iterator j=measured.begin();
        tag_points::const_iterator i=ideal.begin();
  
        fittings::const_iterator bx=basis_x.begin();
        /*fittings::const_iterator by=basis_y.begin();
        fittings::const_iterator bz=basis_z.begin();*/
        std::vector<double> tmp(order*3);
        
        for(;i!=ideal.end();i++, j++, bx++/*, by++, bz++*/)
        {
          
          minc::image3d::IndexType idx;
          error->TransformPhysicalPointToIndex(*i,idx);
          
          tag_point moved;
          for(int k=0;k<3;k++)
          {
            tmp.assign(tmp.size(),0.0);
            for(int t=0;t<order;t++)
              tmp[t+k*order]=(*bx)[t];
            moved[k]=pol_x.fit(tmp,coeff);
          }
          error->SetPixel(idx,(*j).EuclideanDistanceTo(moved));
        }
      }
      save_minc<minc::image3d>(err_f,error);
    }
    
    void load_grid(const char *grid_f,const char *mask_f)
    {
      if(verbose)
        std::cout<<"Loading grid: "<<grid_f<<" and mask: "<<mask_f<<" ... ";
      
      minc::def3d::Pointer grid (minc::def3d::New());
      minc::mask3d::Pointer mask(minc::mask3d::New());
      
      grid=load_minc<minc::def3d>(grid_f);
      mask=load_minc<minc::mask3d>(mask_f);
      
      minc::def3d_iterator it(grid, grid->GetLargestPossibleRegion() );
      int cnt=0;
      for(it.GoToBegin();!it.IsAtEnd();++it)
      {
        tag_point p;
        minc::def3d::IndexType idx=it.GetIndex();
        grid->TransformIndexToPhysicalPoint(idx,p);
        
        minc::mask3d::IndexType idx_m;
        
        if(skip_voxels>1 && ( idx[0]%skip_voxels || idx[1]%skip_voxels || idx[2]%skip_voxels )) continue;
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
    
    void load_pca_matrix(const char * pca_f,int _pcs=-1)
    {
      std::ifstream pca(pca_f);
      pca_matrix.clear();
      
      int ncol=0;
      while(pca.good() && !pca.eof())
      {
        char tmp[65535];
        
        pca.getline(tmp,sizeof(tmp)-1);
        
        if(!strncmp(tmp,"\"PC",3)) //comment
          continue;
        
        if(!strlen(tmp)) break; //eof?
        istringstream ins(tmp);
        std::vector<double> ln;
        char* token;
        for(token=strtok(tmp,",");token!=NULL;token=strtok(NULL,","))
        {
          double val=atof(token);
          ln.push_back(val);
        }
        if(!ncol) ncol=ln.size();
        else if(ncol!=ln.size())
          REPORT_ERROR("Inconsistent number of columns!");
        pca_matrix.push_back(ln);
      }
      
      if(verbose)
        std::cout<<"Loaded "<<pca_matrix.size()<<"x"<<ncol<<" matrix"<<std::endl;
      
      if(_pcs<1) 
        pcs=ncol;
      else 
        pcs=_pcs;
      if(verbose)
        std::cout<<"Going to use "<<pcs<<" PCs for fitting!"<<std::endl;
    }
};

const double Tag_fit::_distance_epsilon=1e-10;

void show_usage (const char * prog)
{
  std::cerr 
      << "Usage: "<<prog<<" <grid1> <mask1> [<grid2> <mask2> .... <grid n> <mask n>] <output.par> " << std::endl
      << "--clobber overwrite files"    << std::endl
      << "--order <n> (3)"<<std::endl
//    << "--keep <part> 0.0-1.0 part of data points to keep (0.8)"<<endl
//    << "--iter <n> maximum number of iterations (200)"<<endl
      << "--remove <pct> 0-1 (0) randomly remove voxels"<<std::endl
      << "--legendre <f> legendre regularization coeefficient"<<std::endl
      << "--error <errr.mnc> output error map"<<std::endl
      << "--pca <matrix.csv> use rotation matrix"<<std::endl
      << "--pcs <N> use this number of PCs"<<std::endl
      << "--cylindrical use cylindrical assumption"<<std::endl
      << "--skip <n> subsample deformation field by factor n"<<endl;
}

int main (int argc, char **argv)
{
  int verbose=0, clobber=0,skip_grid=0,cond=0;
  double max=10.0;
  double keep=1.0;
  double remove=0;
  int iter=1;
  int order=3;
  int skip_voxels=0;
  std::string grid_f,mask_f,output,dump_f;
  std::string residuals_f,error_f,pca_f;
  double legendre=0.0;
  int pcs=-1;
  int cylindrical=0;
  
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"cylindrical", no_argument,   &cylindrical, 1},
    {"order",   required_argument,   0, 'o'},
    {"version", no_argument,         0, 'v'},
    {"legendre",required_argument,   0, 'l'},
    {"error",   required_argument,   0, 'e'},
    {"pca",     required_argument,   0, 'p'},
    {"pcs",     required_argument,   0, 'n'},
    {"skip",    required_argument,   0, 's'},
    {0, 0, 0, 0}
  };
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "o:k:i:vs:l:e:p:n:", long_options, &option_index);

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
      case 'l':
        legendre=atof(optarg);break;
      case 'e':
        error_f=optarg;break;
      case 'p':
        pca_f=optarg;break;
      case 'n':
        pcs=atoi(optarg);break;
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
    Tag_fit fit(basis_functions_x::parameters_no(order),
                CylindricalFunctions::parameters_no(order),
                legendre,max_dev,verbose,cylindrical,skip_voxels);
    
    if(verbose && legendre>0.0)
      std::cout<<"Using legendre coeffecient:"<<legendre<<std::endl;
    
    if(!pca_f.empty() && verbose)
      std::cout<<"Using PCA rotation matrix:"<<pca_f.c_str()<<std::endl;
    
    if(!pca_f.empty())
      fit.load_pca_matrix(pca_f.c_str(),pcs);
    
    output=argv[argc-1];
    
    if (!clobber && !access(output.c_str (), F_OK))
    {
      cerr << output.c_str () << " Exists!" << endl;
      return 1;
    }
    
    for(int i=0;(i+1)<(argc - optind);i+=2)
    {
     
      fit.load_grid(argv[optind+i],argv[optind+i+1]);
    }
    
    
		if(!fit.fit_tags(cond))
		{
			cerr<<"Fitting failed, due to large Standard Deviation !"<<endl;
			return 10;
		}
    if(!error_f.empty())
    {
      fit.save_error(error_f.c_str(),argv[optind]);
    }
    
    std::ofstream cf(output.c_str());
    
    if(!cf.good())
      REPORT_ERROR("Can't open file for writing!");
    
    save_coeff(cf,fit.coeff);
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
