/***************************************************************************
 *   Copyright (C) 2008 by Pierrick COUPE - L. Collins Lab                 *
 *   pcoupe@gonzalo                                                        *
 ***************************************************************************/

/***************************************************************************
 *  This code is the adaptation to minc of the original code               *
 *  by Pierrick COUPE and Pierre YGER                                      *
 *  The original implementation was protected under the licence:           *
 *  IDDN.FR.001.070033.000.S.P.2007.000.21000                              *
 ***************************************************************************/ 

/*                          GENERAL DESCRIPTION                            */
/***************************************************************************
 *  The algorithm for noise removal is described in:                       *
 *                                                                         *
 *  P. Coup�, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
 *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
 *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, * 
 *  Avril 2008                                                             *
 ***************************************************************************/

/*                        ADAPTATION TO SPECKLE                            */
/***************************************************************************
 *  The adaptation to Speckle removal is described in:                     *
 *                                                                         *
 *  P. Coup�, P. Hellier, C. Kervrann, C. Barillot. Bayesian non local     *
 *  means-based speckle filtering. In 5th IEEE International Symposium on  * 
 *  Biomedical Imaging: From Nano to Macro, ISBI'2008, Paris, France,      *
 *  Mai 2008                                                               *
 ***************************************************************************/

/*                      ADAPTATION TO RICAN NOISE                          */
/***************************************************************************
 *  The adaptation to Rician noise is described in:                        *
 *                                                                         *
 *  N. Wiest-Daessl�, S. Prima, P. Coup�, S.P. Morrissey, C. Barillot.     *
 *  Rician noise removal by non-local means filtering for low              *
 *  signal-to-noise ratio MRI: Applications to DT-MRI. In 11th             *
 *  International Conference on Medical Image Computing and                *
 *  Computer-Assisted Intervention, MICCAI'2008, New York, �tats-Unis,     *
 *  Septembre 2008                                                         *
 ***************************************************************************/

#include <iostream>
#include <cstdlib>

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include <minc_1_simple.h> // simple minc reading & writing
#include <time_stamp.h>    // for creating minc style history entry
#include <minc_1_simple_rw.h>
#include "noise_estimate.h"
 
#include <ParseArgv.h>

#include "nl_means_utils.h"
#include "nl_means_block.h"

using namespace std;

int      verbose = FALSE;
int      clobber = FALSE;

double filtering_param   = 0;
double beta    = 1;
int anisotropic  = 0;

double S       = 1;
double M       = 5;

int neighborhoodsize[3]; // size of the patches
int searching[3];        // size of the search area
int testmean   = 1;
int testvar    = 1;
int nb_thread = 4;
double m_min   = 0.95; //threshold for mean test

double v_min   = 0.5; //threshold for var test

int weight_method = 0;
int block      = 1;
int b_space    = 2;
int debug      = FALSE;
int references = FALSE;
char *hallucinate_file = NULL;

static ArgvInfo argTable[] =
{
	{NULL, ARGV_HELP, NULL, NULL,(char*)"###########################################################################"},
	{NULL, ARGV_HELP, NULL, NULL,(char*)"#                                                                         #"},
	{NULL, ARGV_HELP, NULL, NULL,(char*)"#                  Automatic and Multithread Denoising                    #"},
	{NULL, ARGV_HELP, NULL, NULL,(char*)"#                    based on Non Local Means Filter                      #"},
	{NULL, ARGV_HELP, NULL, NULL,(char*)"#                                                                         #"},
	{NULL, ARGV_HELP, NULL, NULL,(char*)"###########################################################################\n"},
    
    
	{NULL, ARGV_HELP, NULL, NULL,(char*)"---------------------------------PARAMETERS--------------------------------"},
	{(char*)"-sigma", ARGV_FLOAT, (char *) 1, (char *) &filtering_param,(char*)"Sigma Value                       [0 = Automatic [default]]"},
	{(char*)"-beta", ARGV_FLOAT, (char *) 1, (char *) &beta,            (char*)"Beta Value                        [default 1]"},
	{(char*)"-v", ARGV_FLOAT, (char *) 1, (char *) &S,                  (char*)"Neighboring size : \n\t\t 1 : 26 neighbors [default] \n\t\t 2 : 124 neighbors ..."},
	{(char*)"-d", ARGV_FLOAT, (char *) 1, (char *) &M,                  (char*)"Search Volume size                [default 5 : 1331 neighbors]"},
	{(char*)"-w", ARGV_INT, (char *) 1, (char *) &weight_method,        (char*)"Weighting method : \n\t\t 0 : L2-norm (Gaussian noise) [default]\n\t\t 1 : Pearson Divergence (Speckle)\n\t\t 2 : L2-norm + Bais correction (Rician noise)"},
	{(char*)"-aniso", ARGV_CONSTANT, (char *) 1, (char *) &anisotropic, (char*)"The neighborhood size is adapted on the anisotropy of the image, else the neighborhood is cubic \n"},
	{NULL, ARGV_HELP, NULL, NULL,(char*)"----------------------------------OPTIONS----------------------------------"},
	{(char*)"-block", ARGV_INT, (char *) 1, (char *) &block,       (char*)"NL-means based on block approach  [default 1]"},
	{(char*)"-b_space", ARGV_INT, (char *) 1, (char *) &b_space,   (char*)"Distance between blocks           [default 2]"},
	{(char*)"-m_min", ARGV_FLOAT, (char *) 1, (char *) &m_min,     (char*)"Lowest bound of mean ratio        [default 0.95]"},
	{(char*)"-v_min", ARGV_FLOAT, (char *) 1, (char *) &v_min,     (char*)"Lowest bound of variance ratio    [default 0.5]"},
	{(char*)"-mt", ARGV_INT, (char *) 1, (char *) &nb_thread,      (char*)"Number of thread                  [default 4]"},
  {(char*)"-hallucinate",  ARGV_STRING,  (char*)1,  (char*) &hallucinate_file,     (char*)"Hallucinate this file (experimental)."},
    
	{NULL, ARGV_HELP, NULL, NULL,(char*)"\n---------------------------------------------------------------------------"},
	{(char*)"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,(char*)"Print out extra information."},
	{(char*)"-debug", ARGV_CONSTANT, (char *)TRUE, (char *)&debug,	(char*)"Print out  even more information."},
	{(char*)"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,(char*)"Clobber existing files.\n"},
  {(char*)"-references",ARGV_CONSTANT, (char *)TRUE, (char *)&references,(char*)"Print citation references\n"},
    
	{NULL, ARGV_END, NULL, NULL, (char*)"       "},
	{NULL, ARGV_END, NULL, NULL, (char*)"       "},
	{NULL, ARGV_END, NULL, NULL, NULL}
};

void print_references(void) 
{
 std::cout<<"\
 \
\
                          GENERAL DESCRIPTION                            \n\
 \n\
   The algorithm for noise removal is described in:                       \n\
                                                                          \n\
   P. Coup�, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     \n\
   An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic\n\
   Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441,  \n\
   Avril 2008                                                             \n\
 \n\
\n\
                        ADAPTATION TO SPECKLE                            \n\
\n\
   The adaptation to Speckle removal is described in:                     \n\
                                                                          \n\
   P. Coup�, P. Hellier, C. Kervrann, C. Barillot. Bayesian non local     \n\
   means-based speckle filtering. In 5th IEEE International Symposium on   \n\
   Biomedical Imaging: From Nano to Macro, ISBI'2008, Paris, France,      \n\
   Mai 2008                                                               \n\
 \n\
\n\
                      ADAPTATION TO RICAN NOISE                          \n\
\n\
   The adaptation to Rician noise is described in:                        \n\
                                                                          \n\
   N. Wiest-Daessl�, S. Prima, P. Coup�, S.P. Morrissey, C. Barillot.     \n\
   Rician noise removal by non-local means filtering for low              \n\
   signal-to-noise ratio MRI: Applications to DT-MRI. In 11th             \n\
   International Conference on Medical Image Computing and                \n\
   Computer-Assisted Intervention, MICCAI'2008, New York, �tats-Unis,     \n\
   Septembre 2008                                                         \n\
";
}

// Function to perform preprocessing and call the denoising function
void Exec(minc::simple_volume<float> & in, float *ima_out, int* vol_size, Real * vol_res, float *hallucinate_in)
{
  float *ima_in=in.c_buf();

	if(vol_size[2] < 2*nb_thread)
	{
		std::cout << "\n------------------------------------------------" << std::endl;
		std::cout << "!The number of slices is too small (<nb_thread)!" << std::endl;
		std::cout << "!!    => Set voxelwise mode (block = 0)       !!" << std::endl;
		std::cout << "------------------------------------------------"<< std::endl;
		block =0;
	}

	if (anisotropic == 0)
	{
		neighborhoodsize[0] = (int) S;
		neighborhoodsize[1] = (int) S;
		neighborhoodsize[2] = (int) S;
		searching[0] = (int) M;
		searching[1] = (int) M;
		searching[2] = (int) M;
	}
	else
	{
    if(verbose) {
      std::cout << "\n------------------------------------------------" << std::endl;
      std::cout << "                 Anisotropy mode                 " << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << " The neighborhood size is adapted to the anisotropy of the image" << std::endl;
    }
    
    // set default parameters
    neighborhoodsize[0] = (int) S;  //2
    neighborhoodsize[1] = (int) S;  //2
    neighborhoodsize[2] = (int) S;  //2
    searching[0] = (int) M;    //7
    searching[1] = (int) M;    //7
    searching[2] = (int) M;   //7

    std::cout<<" -- Neighborhoodsize: "<<neighborhoodsize[0]<<","<<neighborhoodsize[1]<<","<<neighborhoodsize[2]<<std::endl;
    std::cout<<" -- Searching: "<<searching[0]<<","<<searching[1]<<","<<searching[2]<<std::endl;
    
    //check 2D dimension
    if (vol_res[0]>vol_res[1] && vol_res[0]>vol_res[2])
    {
      neighborhoodsize[0] = (int) 0;
      searching[0] = (int) 1;
    }
    else if (vol_res[1]>vol_res[2] && vol_res[1]>vol_res[0])
    {
      neighborhoodsize[1] = (int) 0;
      searching[1] = (int) 1;
    }
    else if (vol_res[2]>vol_res[1] && vol_res[2]>vol_res[0])
    {
      neighborhoodsize[2] = (int)0;
      searching[2] = (int) 1;
    }

    std::cout<<" -- Resolution: "<<vol_res[0]<<","<<vol_res[1]<<","<<vol_res[2]<<std::endl;
    std::cout<<" -- Neighborhoodsize: "<<neighborhoodsize[0]<<","<<neighborhoodsize[1]<<","<<neighborhoodsize[2]<<std::endl;
    std::cout<<" -- Searching: "<<searching[0]<<","<<searching[1]<<","<<searching[2]<<std::endl;
  }
  
  if(verbose) {
    std::cout <<"\n--------------------------------------------------"<< std::endl;
    std::cout <<    "                   Parameters                     " << std::endl;
    std::cout <<"--------------------------------------------------\n"<< std::endl;
   
    std::cout <<      "Sigma           : ";
    if (filtering_param == 0)
      std::cout << "Automatic: only Rician or Gaussian noise" << std::endl;
    else
      std::cout << filtering_param << std::endl;
    if (beta!=1.0)
      std::cout << "Weighting parameter " << beta << std::endl;
      
    std::cout <<      "Ni              : " << 2*neighborhoodsize[0]+1 << " x " << 2*neighborhoodsize[1]+1 << " x " << 2*neighborhoodsize[2]+1 << std::endl;
    std::cout <<      "Vi              : " << 2*searching[0]+1 << " x " << 2*searching[1]+1 << " x " << 2*searching[2]+1 << std::endl;
      
    if (testmean)
      std::cout <<    "Mean Test       : Yes, " << m_min << " < X < " << 1/m_min << std::endl;
    else
      std::cout << "Mean Test       : No" << std::endl;
      
    if (testvar)
      std::cout <<    "Variance Test   : Yes, " << v_min << " < X < " << 1/ v_min << std::endl;
    else
      std::cout << "Variance Test   : No" << std::endl;
      
    std::cout << std::endl;
      
    if (anisotropic == 1)
      std::cout <<  "Anisotropic mode activated       : Yes" << std::endl;
    else
      std::cout <<  "Anisotropic mode activated       : No" << std::endl;
  } 
  
	if (block == 1)
	{
    if(verbose) {
      std::cout << "Block Implementaiton of NL-means : Yes" << std::endl;
      std::cout << "--> Distance between blocks      : " << b_space << std::endl;
    }
		if (neighborhoodsize[0] <(b_space/2))
		{
			std::cout << "We must have Ni < (D_block)/2" << std::endl;
			exit(0);
		}
    
		if (anisotropic == 1)
		{
			std::cout << "Can't activate isotropic mode in block version" << std::endl;
			exit(0);
		}
    
	}
  
  if(verbose) {
    if (block == 0)
      std::cout << "Block Implementation of NL-means : No " << std::endl;
    if (weight_method == 0)
      std::cout << "Weighting Method                 : L2-norm (Gaussian noise) " << std::endl;
    if (weight_method == 1)
      std::cout << "Weighting Method                 : Pearson Divergence (Speckle)" << std::endl;
    if (weight_method == 2)
      std::cout << "Weighting Method                 : L2-norm + Bais correction (Rician noise) " << std::endl;
    
    std::cout<< "\n" << std::endl;
  }
  
	/* local and local variance storage */
	float *mean_map, *var_map;
	mean_map = new float[vol_size[0]*vol_size[1]*vol_size[2]];
	var_map = new float[vol_size[0]*vol_size[1]*vol_size[2]];
	for (int i = 0; i < vol_size[0] *  vol_size[1] * vol_size[2]; i++)
	{
		mean_map[i] = 0.0;
		var_map[i] = 0.0;
	}

  if(verbose) {
    std::cout <<"\n--------------------------------------------------"<< std::endl;
    std::cout <<    "                  Preprocessing                      " << std::endl;
    std::cout <<"--------------------------------------------------\n"<< std::endl;
  }
  
   // Computation of the local means
	Preprocessing(ima_in,mean_map,neighborhoodsize,vol_size);
    
	if (testvar == 1)
      // Computation of the local variances
		Preprocessing2(ima_in,mean_map,var_map,neighborhoodsize,vol_size);
  
	if (filtering_param == 0)
	{
    
    if(weight_method ==1)
    {
      std::cout <<" - ERROR: automatic variance estimation is not available for Spekle Noise:  option '-w 1'"<<std::endl
                <<"          Use option -sigma"<<std::endl;
      return;
    }
    double mean_val=0.0;
    filtering_param=minc::noise_estimate(in,mean_val,weight_method==0,verbose);//use gaussian estimate appropriately
	} 

  if(verbose) {  
    std::cout <<"\n--------------------------------------------------"<< std::endl;
    std::cout <<    "                  Denoising                      " << std::endl;
    std::cout <<"--------------------------------------------------\n"<< std::endl;
  }
  
	if (block==0)   
		denoise_mt(ima_in,ima_out,mean_map,var_map, filtering_param,beta,
        neighborhoodsize,searching,
        testmean,testvar,m_min,v_min,weight_method,vol_size,hallucinate_in);
  
	if (block==1)  
		denoise_block_mt(ima_in,ima_out,mean_map,var_map,filtering_param,beta,
        neighborhoodsize,searching,
        testmean,testvar,m_min,v_min,weight_method,b_space,vol_size,hallucinate_in);

	delete [] mean_map;
	delete [] var_map;
}


// Main function to read/write the volume
int main(int argc, char *argv[])
{

  //Volume  in_vol, in_volio_float, out_volio_float; // All the work is perfomed in float
	float   *in_vol_float, *out_vol_float,*in_hallucinate=NULL;
  minc::simple_volume<float> in_vol;

	char    *in_file;
	char    *out_file;
	int      c, in_ndims;
	float     min, max;
	Status status = ERROR;
	nc_type  datatype;
	BOOLEAN  signed_flag;
	minc_input_options in_ops;
	Real     min_value, max_value;
	int vol_size[3] = {0, 0, 0}; /* size */
	Real vol_res[3] = {0.0, 0.0, 0.0};  /* resolutions */
  
	char *history = time_stamp(argc, argv); //maybe we should free it afterwards

  std::cout<<"\
###########################################################################\n\
#                                                                         #\n\
#                  Automatic and Multithread Denoising                    #\n\
#                    based on Non Local Means Filter                      #\n\
#                                                                         #\n\
###########################################################################\n\
\n\
Copyright (C) 2008 by Pierrick COUPE - L. Collins Lab \n\
\n\
This code is the adaptation to minc of the original code \n\
by Pierrick COUPE and Pierre YGER\n\
The original implementation was protected under the licence: \n\
IDDN.FR.001.070033.000.S.P.2007.000.21000\n\n";
  
	/* get args */
	if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2))
	{
		fprintf(stderr,"\nUsage: %s [<options>] <infile.mnc> <outfile.mnc>\n", argv[0]);
		fprintf(stderr, "       %s [-help]\n\n", argv[0]);
    if(references) print_references();
		exit(EXIT_FAILURE);
	}
  
  if(references) print_references();

	in_file = argv[1];
	out_file = argv[2];
  
	/* check for the infile and outfile */
  /*if(access(in_file, F_OK) != 0){
	fprintf(stderr, "%s: Couldn't find %s\n\n", argv[0], in_file);
	exit(EXIT_FAILURE);
}*/
	if(access(out_file, F_OK) == 0 && !clobber){
		fprintf(stderr, "%s: %s exists! (use -clobber to overwrite)\n\n", argv[0], out_file);
		exit(EXIT_FAILURE);
	}
	try {
      
		minc::minc_1_reader minc_reader;
		minc_reader.open(in_file);
     
		/* read in the input file */
     //in_ndims = get_minc_file_n_dimensions(in_file);
     //set_default_minc_input_options(&in_ops);


		if((in_ndims=minc_reader.dim_no()) != 3){
			fprintf(stderr, "%s: Only 3D volume is supported for the moment\n", argv[0]);
			return EXIT_FAILURE;
		}
     
		vol_res[0] = ABS(minc_reader.nspacing(1));
		vol_res[1] = ABS(minc_reader.nspacing(2));
		vol_res[2] = ABS(minc_reader.nspacing(3));
		vol_size[0] = minc_reader.ndim(1);
		vol_size[1] = minc_reader.ndim(2);
		vol_size[2] = minc_reader.ndim(3);
    
    
		//in_vol_float =  new float[vol_size[0]*vol_size[1]*vol_size[2]];
    
		out_vol_float = new float[vol_size[0]*vol_size[1]*vol_size[2]];
    
    /*
		minc_reader.setup_read_float();
		minc::load_standard_volume<float>(minc_reader,in_vol_float);*/
    minc::load_simple_volume(minc_reader,in_vol);
    in_vol_float=in_vol.c_buf();
    
    if(hallucinate_file)
    {
      minc::minc_1_reader h_reader;
      h_reader.open(hallucinate_file);
      for(int i=1;i<4;i++)
        if(h_reader.ndim(i)!=minc_reader.ndim(i))
        {
          fprintf(stderr, "%s: Inconsistent dimension size of hallucination file\n", argv[0]);
          return EXIT_FAILURE;
        }
      in_hallucinate=new float[vol_size[0]*vol_size[1]*vol_size[2]];
      h_reader.setup_read_float();
      minc::load_standard_volume<float>(h_reader,in_hallucinate);
    }
    
		min_value=1e10;
		max_value=-1e10;
    //calculate min & max values
    
		for (int i = 0; i < vol_size[0] *  vol_size[1] * vol_size[2]; i++)
		{
			if(min_value>in_vol_float[i]) min_value=in_vol_float[i];
			if(max_value<in_vol_float[i]) max_value=in_vol_float[i];
			out_vol_float[i] = 0.0;
		}
    
    /*if(filtering_param==0.0) //estimating noise
    {
      filtering_param=noise_estimate(in_vol,mean_value,false,verbose);
    }*/
    
		if(verbose)
		{
      
			fprintf(stdout, " | Input file:     %s\n", in_file);
			fprintf(stdout, " | Input ndims:    %d\n", in_ndims);
			fprintf(stdout, " | min/max:        [%8.3f:%8.3f]\n", min_value, max_value);
      
      if(hallucinate_file) 
        fprintf(stdout, " | Hallucination file: %s\n", hallucinate_file);
			fprintf(stdout, " | Output files:\n");
		}
      
		if (weight_method == 1) 
		{
			if (filtering_param == 0.0)
			{
				fprintf(stderr, "You have to give a sigma value\n\n");
				return EXIT_FAILURE;
			}
		}
    
    if(verbose)
    {
      std::cout <<"\n--------------------------------------------------"<< std::endl;
      std::cout <<  "|                                                |" << std::endl;
      std::cout <<  "|       3D NonLocal means denoising filter       |" << std::endl;
      std::cout <<  "|                                                |" << std::endl;
      std::cout <<  "--------------------------------------------------" << std::endl;
    }
       
		if (weight_method != 2)
		{
      /*min = FindMinVolume(in_vol_float, vol_size);
      max = FindMaxVolume(in_vol_float, vol_size);*/
      // Linear scaling of intensities between 0 -255 to keep homogenous parameters 
      //LinearScaling(in_vol_float, min, max, vol_size);
      
      if(in_hallucinate) //scale hallucination volume too
      {
        /*min = FindMinVolume(in_hallucinate, vol_size);
        max = FindMaxVolume(in_hallucinate, vol_size);*/
        // Linear scaling of intensities between 0 -255 to keep homogenous parameters 
        //LinearScaling(in_hallucinate, min, max, vol_size);
      }
      
      /*if(verbose)
      {
        std::cout << "\n--------------------------------------------------" << std::endl;
        std::cout << "Linear scaling of the intensities between 0 - 255 " << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
      }*/
      
      ///@todo rescale the sigma of noise!!!
		}
    
    
    // Excecusion of the program
		Exec(in_vol, out_vol_float, vol_size, vol_res, in_hallucinate);

		if (weight_method != 2)
		{
      // Linear resclaing of intensities between min and max
			//LinearReScaling(out_vol_float, min, max, vol_size);
      /*
      if(verbose) {
        std::cout << "\n--------------------------------------------------" << std::endl;
        std::cout << "    Linear rescaling of the intensities between "<<min<<" - "<<max<< std::endl;
        std::cout << "-------------------------------------------------\n" << std::endl;
      }*/
		}
    


		minc::minc_1_writer minc_writer;
		minc_writer.open(out_file,minc_reader.info(),2,minc_reader.datatype());
		minc_writer.copy_headers(minc_reader);
		minc_writer.append_history(history);
      
		minc_writer.setup_write_float();
		minc::save_standard_volume<float>(minc_writer,out_vol_float);
      
        
		//delete [] in_vol_float;
		delete [] out_vol_float;
		
	} catch (const minc::generic_error & err) {
		std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
		std::cerr << err.msg()<<std::endl;
		return 1;
	}
	return (EXIT_SUCCESS);
 
}

