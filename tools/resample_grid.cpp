#include "minc_wrappers.h"
#include <minc_helpers.h>

#include <unistd.h>
#include <getopt.h>
//#include "data_proc.h"

//#include <itkVectorLinearInterpolateImageFunction.h>
#include "mincVectorBSplineInterpolate.h"
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_det.h>
#include <vnl/vnl_vector_fixed_ref.h>
#include <vnl/vnl_vector.h>

typedef minc::mincVectorBSplineInterpolate<minc::def3d,double> Vector_Interpolator;

#include <volume_io.h>

using namespace std;
using namespace minc;

void show_usage (const char * prog)
{
  std::cerr 
    << "Usage: "<<prog<<" <input> <xfm> <output.mnc> " << std::endl
    << "--clobber overwrite files"    << std::endl
    << "--like <example> (default behaviour analogous to use_input_sampling)"<<std::endl
    << "--order <n> spline order 0-5 ; 1st - linear interpolation, default 3"<<std::endl;
}

int main (int argc, char **argv)
{
  int verbose=0, clobber=0,skip_grid=0,order=3;
  std::string like_f,xfm_f,output_f,input_f;
  
  static struct option long_options[] = {
		{"verbose", no_argument,       &verbose, 1},
		{"quiet",   no_argument,       &verbose, 0},
		{"clobber", no_argument,       &clobber, 1},
		{"like",    required_argument, 0, 'l'},
    {"order",    required_argument, 0, 'o'},
		{0, 0, 0, 0}
		};
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "vl:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
			{
			case 0:
				break;
			case 'v':
				cout << "Version: 0.1" << endl;
				return 0;
      case 'l':
        like_f=optarg; break;
      case 'o':
        order=atoi(optarg);
        break;
			case '?':
				/* getopt_long already printed an error message. */
			default:
				show_usage (argv[0]);
				return 1;
			}
    }

	if ((argc - optind) < 3) {
		show_usage(argv[0]);
		return 1;
	}
  input_f=argv[optind];
  xfm_f=argv[optind+1];
  output_f=argv[optind+2];
  
  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }
	try
  {
    const double delta=1e-5;
		minc::def3d::Pointer in_grid(minc::def3d::New());
    minc::def3d::Pointer out_grid(minc::def3d::New());
    load_minc(input_f.c_str(),in_grid);
    
    if(!like_f.empty())
      imitate_minc(like_f.c_str(),out_grid);
    else
      allocate_same(out_grid,in_grid);
    
    Vector_Interpolator::Pointer interpolator(Vector_Interpolator::New());
    interpolator->SetSplineOrder(order);
    
    interpolator->SetInputImage(in_grid);
    
    General_transform xfm;
   /* read in  the input transformation */
    if(input_transform_file((char*)xfm_f.c_str(), &xfm)){
      std::cerr<<"Error reading in xfm"<<xfm_f.c_str()<< std::endl;
      return 1;
    }
    tag_point zero_vec;
    zero_vec[0]=zero_vec[1]=zero_vec[2]=0;
    
    minc::def3d_iterator it(out_grid,out_grid->GetLargestPossibleRegion());
    vnl_vector<double> _in(3);
    for(it.GoToBegin(); !it.IsAtEnd();++it)
    {
      itk::ContinuousIndex< double, 3 > cidx; 
      tag_point dst,orig;
      Vector_Interpolator::OutputType in;
      vnl_vector_fixed<double,3> _out;
      double u1,v1,w1;
      double u2,v2,w2;
      vnl_matrix_fixed< double, 3, 3 > m;
      out_grid->TransformIndexToPhysicalPoint(it.GetIndex(),dst);
      general_inverse_transform_point(&xfm,dst[0],dst[1],dst[2],&u1,&v1,&w1);
      orig[0]=u1;orig[1]=v1;orig[2]=w1;
      if(!in_grid->TransformPhysicalPointToContinuousIndex<double>(orig,cidx))
      {
        it.Value()[0]=0;
        it.Value()[1]=0;
        it.Value()[2]=0;
        continue;
      }
      _in.set(interpolator->Evaluate(orig).GetDataPointer());
      //estimate 1st derivatives
      general_transform_point(&xfm, orig[0]-delta, orig[1], orig[2],&u1, &v1, &w1);
      general_transform_point(&xfm, orig[0]+delta, orig[1], orig[2],&u2, &v2, &w2);
      m(0,0)=(u2-u1)/(2*delta);
      m(0,1)=(v2-v1)/(2*delta);
      m(0,2)=(w2-w1)/(2*delta);
      
      general_transform_point(&xfm, orig[0], orig[1]-delta, orig[2],&u1, &v1, &w1);
      general_transform_point(&xfm, orig[0], orig[1]+delta, orig[2],&u2, &v2, &w2);
      m(1,0)=(u2-u1)/(2*delta);
      m(1,1)=(v2-v1)/(2*delta);
      m(1,2)=(w2-w1)/(2*delta);
      
      general_transform_point(&xfm, orig[0], orig[1], orig[2]-delta,&u1, &v1, &w1);
      general_transform_point(&xfm, orig[0], orig[1], orig[2]+delta,&u2, &v2, &w2);
      m(2,0)=(u2-u1)/(2*delta);
      m(2,1)=(v2-v1)/(2*delta);
      m(2,2)=(w2-w1)/(2*delta);
      
      //_out=m * _in ;
      _out=_in*m  ;
      //it.Value()=	vnl_det(m)-1;
      //vnl_vector_fixed_ref<float,3) _in(cnv.
      it.Value()[0]=_out[0];
      it.Value()[1]=_out[1];
      it.Value()[2]=_out[2];
      
    }
    delete_general_transform(&xfm);

    save_minc<def3d>(output_f.c_str(),out_grid);
    
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
