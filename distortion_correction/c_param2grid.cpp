/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
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
#include <iostream>
#include <fstream>
#include "mincMeanSquaresImageToImageMetric.h"
#include "sphericalHarmonicsTransform.h"

#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>

using namespace std;
using namespace minc;

void show_usage (void)
{
  std::cerr 
    << "Usage: param2grid <parameters> <output_grid.mnc> " << endl
    << "--clobber overwrite files"    << endl
    << "--spacing <n> spacing of grid file in mm, default 4"<< endl
    << "--max <dist> maximum distance, default 5mm"<<endl
    << "--extent <mm> - data span ( -extent/2 : extent/2 in all dimensions), default 300mm"<<endl;

}


int main (int argc, char **argv)
{
  int verbose=0, clobber=0,skip_grid=0;
  double max=5.0;
  double extent=300;

  static struct option long_options[] = {
		{"verbose", no_argument,       &verbose, 1},
		{"quiet",   no_argument,       &verbose, 0},
		{"clobber", no_argument,       &clobber, 1},
		{"spacing", required_argument, 0, 's'},
		{"max",     required_argument, 0, 'm'},
		{"version", no_argument,       0, 'v'},
    {"extent", required_argument,  0, 'e'},
		{0, 0, 0, 0}
		};
  
  double spacing=4.0;
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "s:m:v", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
			{
			case 0:
				break;
			case 's':
				spacing=atof(optarg);
				break;
			case 'v':
				cout << "Version: 1.0" << endl;
				return 0;
      case 'm':
        max=atof(optarg); break;
      case 'e':
        extent=atof(optarg); break;
			case '?':
				/* getopt_long already printed an error message. */
			default:
				show_usage ();
				return 1;
			}
    }

	if ((argc - optind) < 2) {
		show_usage ();
		return 1;
	}
  std::string input=argv[optind];
  std::string output=argv[optind+1];
	try
  {
    gsl_rng_env_setup();
    
		typedef minc::CylindricalHarmonicsTransform TransformType;
		TransformType::ParametersType finalParameters;
		load_parameters(input.c_str(),finalParameters);
		TransformType::Pointer finalTransform = TransformType::New();
    cout<<"Loaded parameters:"<<finalParameters<<endl;
		finalTransform->ImportParameters( finalParameters , true);
    
    cout<<"Imported!"<<endl;
		minc::def3d::Pointer grid(minc::def3d::New());
		allocate_image3d(grid, fixed_vec<3, unsigned int>(extent/spacing), fixed_vec<3, double>(spacing), fixed_vec<3, double>(-extent/2));
		
		if(verbose) 
		{
			cout<<"Generating a grid file, ";
			cout<<"extent: "<<extent<<" spacing: "<<spacing<<" ..."<<std::flush;
		}
		
    def3d_iterator it(grid,grid->GetLargestPossibleRegion());
		for(it.GoToBegin();!it.IsAtEnd();++it) {
      tag_point p,p2;
      grid->TransformIndexToPhysicalPoint(it.GetIndex(),p);
			p2=finalTransform->TransformPointUnCached(p);
      def_vector moved;
			moved[0]=p2[0]-p[0];
			moved[1]=p2[1]-p[1];
			moved[2]=p2[2]-p[2];
      if(fabs(moved[0])>max || fabs(moved[1])>max ||fabs(moved[2])>max)
        moved[0]=moved[1]=moved[2]=0.0;
			
      it.Value()=moved;
    }
		if(verbose)
			cout<<"Done!"<<endl;
    save_minc<minc::def3d>(output.c_str(), grid);
		
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
