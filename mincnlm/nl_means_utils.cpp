#include "nl_means_utils.h"

extern int      verbose;
extern int      debug;

float FindMinVolume(float *ima_in, int *vol_size)
{

	float min, tmp;
 
	min = ima_in[0];


	for (int i = 1; i < vol_size[0] * vol_size[1] * vol_size[2]; i++)
	{
       // tmp = get_volume_real_value(ima_in, i, j, k, 0, 0)
		tmp = ima_in[i];
		if (tmp < min) min = tmp;
	}
   
  
	return min;
}

float FindMaxVolume(float *ima_in, int *vol_size)
{

	float max, tmp;
 
	max = ima_in[0];;

	for (int i = 1; i < vol_size[0] * vol_size[1] * vol_size[2]; i++)
	{
        
		tmp = ima_in[i];
		if (tmp > max) max = tmp;
	}
  
	return max;
}

// Linear scaling between 0 - 255
void LinearScaling(float *ima_in, double mini, double maxi, int *vol_size)
{

	double a = 255/(maxi-mini);
	double b = 255*mini/(maxi-mini);

	for (int i = 0; i < vol_size[0]  * vol_size[1] * vol_size[2]; i++)
	{
        //double tmp = get_volume_real_value(ima_in, i, j, k, 0, 0)*a-b;
        //set_volume_real_value(ima_in,i,j,k,0,0, (float)tmp);
		double tmp = ima_in[i]*a-b;
		ima_in[i] = (float)tmp;
	}

}

// Linear rescaling between mini - maxi
void LinearReScaling(float *ima_in, double mini, double maxi, int *vol_size)
{
	for (int i = 0; i < vol_size[0] * vol_size[1] * vol_size[2]; i++)
	{
        //double tmp = get_volume_real_value(ima_in, i, j, k, 0, 0) *(maxi-mini)/255 + mini;
        //set_volume_real_value(ima_in,i,j,k,0,0, std::min(maxi,std::max(mini,tmp)) );
        
		double tmp = ima_in[i] *(maxi-mini)/255 + mini;
		ima_in[i] = std::min(maxi,std::max(mini,tmp));
        
	}
  
}


float Neighborhood_Mean(float *ima_in,int x,int y,int z,int *neighborhoodsize, int * vol_size)
{
	int x_pos,y_pos,z_pos;
	bool is_outside;
	float global_sum = 0;
	int nb_inside = 0;
 
  
	for (int c = 0; c<(2*neighborhoodsize[2]+1);c++){
		for (int b = 0; b<(2*neighborhoodsize[1]+1);b++){
			for (int a = 0;a<(2*neighborhoodsize[0]+1);a++){

				is_outside = false;
				x_pos = x+a-neighborhoodsize[0];
				y_pos = y+b-neighborhoodsize[1];
				z_pos = z+c-neighborhoodsize[2];

				if ((z_pos < 0) || (z_pos > vol_size[2]-1)) is_outside = true;
				if ((y_pos < 0) || (y_pos > vol_size[1]-1)) is_outside = true;
				if ((x_pos < 0) || (x_pos > vol_size[0]-1)) is_outside = true;
				if (!is_outside)
				{
          //global_sum = global_sum + get_volume_real_value(ima_in,x_pos,y_pos,z_pos,0,0);
					global_sum = global_sum + ima_in[z_pos*(vol_size[0]*vol_size[1])+(y_pos*vol_size[0])+x_pos];
					nb_inside++;
				}
			}
		}
	}
	float mean = global_sum/nb_inside;
	return mean;
}

void Preprocessing(float *ima_in, float *mean_map,int *neighborhoodsize, int *vol_size)
{

  if(verbose)
    std::cout << "Preprocessing of the volume (mean)..." << std::endl;

  
	for (int k = 0; k < vol_size[2]; k++) {

		for (int j = 0; j < vol_size[1];j++) {

			for (int i = 0; i < vol_size[0];i++) {
        
        //set_volume_real_value(mean_map,i,j,k,0,0,Neighborhood_Mean(ima_in,i,j,k,neighborhoodsize, vol_size));
				mean_map[k*(vol_size[0]*vol_size[1])+(j*vol_size[0])+i] = Neighborhood_Mean(ima_in,i,j,k,neighborhoodsize, vol_size);
        
			}
		}
	}
}



float Neighborhood_Var(float *ima_in, float *ima_mean, int x, int y, int z,int *neighborhoodsize, int * vol_size)
{
  
	int x_pos,y_pos,z_pos;
	bool is_outside;

	float value =0.0;
	float global_sum = 0.0;
	int nb_inside = 0;
  
	int POS = 0;

	for (int c = 0; c<(2*neighborhoodsize[2]+1);c++){
    
		for (int b = 0; b<(2*neighborhoodsize[1]+1);b++){
      
			for (int a = 0;a<(2*neighborhoodsize[0]+1);a++){

				is_outside = false;
				value = 0.0;
				x_pos = x+a-neighborhoodsize[0];
				y_pos = y+b-neighborhoodsize[1];
				z_pos = z+c-neighborhoodsize[2];

				if ((z_pos < 0) || (z_pos > vol_size[2]-1)) is_outside = true;
				if ((y_pos < 0) || (y_pos > vol_size[1]-1)) is_outside = true;
				if ((x_pos < 0) || (x_pos > vol_size[0]-1)) is_outside = true;
				if (!is_outside) {
         // value = get_volume_real_value(ima_in,x_pos,y_pos,z_pos,0,0) - get_volume_real_value(ima_mean,x_pos,y_pos,z_pos,0,0);
					value =  ima_in[z_pos*(vol_size[0]*vol_size[1])+(y_pos*vol_size[0])+x_pos] -  ima_mean[z_pos*(vol_size[0]*vol_size[1])+(y_pos*vol_size[0])+x_pos];
					value = value * value;
					global_sum = global_sum + value ; 
					nb_inside++;
				}
			}
		}
	}
	float var = global_sum/(nb_inside-1);
	return var;
}


void Preprocessing2(float *ima_in, float *mean_map, float *var_map,int *neighborhoodsize, int * vol_size)
{
  if(verbose)
    std::cout << "Preprocessing of the volume (variance)..." << std::endl;

  
	for (int k = 0; k < vol_size[2]; k++) 
	{
		for (int j = 0; j < vol_size[1];j++) 
		{
			for (int i = 0; i < vol_size[0];i++) 
			{
        //set_volume_real_value(var_map,i,j,k,0,0,Neighborhood_Var(ima_in,mean_map,i,j,k,neighborhoodsize, vol_size));
				var_map[k*(vol_size[0]*vol_size[1])+(j*vol_size[0])+i] = Neighborhood_Var(ima_in,mean_map,i,j,k,neighborhoodsize, vol_size);
			}
		}
	}
}

//STd estimation of the gaussian noise by pseudo residu
void Variance_Estimation_1(float *ima_in, int weight_method, double &s_constant, int * vol_size)
{
 
	float value = 0.0;
	float global_variance = 0.0;
	float global_mean     = 0.0;

 
	for (int k = 0; k < vol_size[2]; k++) {
   
		for (int j = 0; j < vol_size[1];j++) {
    
			for (int i = 0; i < vol_size[0];i++) {
        
        
    
				value = 6*ima_in[k*(vol_size[0]*vol_size[1])+(j*vol_size[0])+i] - (
              
        
						ima_in[k*(vol_size[0]*vol_size[1])+(j*vol_size[0])+std::max(0,i-1)] +
            
          
						ima_in[k*(vol_size[0]*vol_size[1])+(j*vol_size[0])+ std::min(vol_size[0]-1,i+1)] +
            
           
						ima_in[k*(vol_size[0]*vol_size[1])+(std::max(0,j-1)*vol_size[0])+ i] +
              
           
						ima_in[k*(vol_size[0]*vol_size[1])+(std::min(vol_size[1]-1,j+1)*vol_size[0])+ i] +
            
           
						ima_in[std::max(0,k-1)*(vol_size[0]*vol_size[1])+(j*vol_size[0])+ i] +
            
           
						ima_in[std::min(vol_size[2]-1,k+1)*(vol_size[0]*vol_size[1])+(j*vol_size[0])+ i] );

        
				global_variance = global_variance + (value*value)/42;
      
				global_mean     = global_mean + ima_in[k*(vol_size[0]*vol_size[1])+(j*vol_size[0])+ i];
			}
		}
	}

	global_variance = global_variance / (vol_size[0]*vol_size[1]*vol_size[2]);
	global_mean     = global_mean     / (vol_size[0]*vol_size[1]*vol_size[2]);


// if ((weight_method == 0) ||(weight_method == 2))
//   {
//     double  STDcorrection = 1.; //Bais correction of the std estimation for Gaussian denoising
	//                           
	// 
//     STDcorrection = 0.0417 * sqrt(global_variance) + 0.0208; //empirical correction based on brainweb experiments between 0 - 255
	// 
//     if ( STDcorrection <1 )
//     {
//       s_constant = sqrt(2*global_variance*s_constant*STDcorrection);
//       //std::cout<<"Correction factor: "<<STDcorrection<<std::endl;
//     }
//   else  s_constant = sqrt(2*global_variance*s_constant);
	// 
//     std::cout << "\nEstimation of the noise standard deviation: " << sqrt(global_variance)  <<std::endl;
//     std::cout << "--> Smoothing paramter h = " << s_constant << std::endl;
//   }

	if (block == 0)
		s_constant = sqrt(global_variance*s_constant);
	else 
		s_constant = sqrt(0.5*global_variance*s_constant);

  if(verbose) {
    std::cout << "\nEstimation of the noise standard deviation: " << sqrt(global_variance)  <<std::endl;
    std::cout << "--> Smoothing paramter h = " << s_constant << std::endl;
  }
}
