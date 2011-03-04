#include "nl_means_fft.h"

extern int      verbose;
extern int      debug;

//Computation of the SS of vector
float L2_norm(float *V1, int *neighborhoodsize)
{

  float result;
  result = 0.0;
  int size = (2*neighborhoodsize[0]+1)*(2*neighborhoodsize[1]+1)*(2*neighborhoodsize[2]+1);
  for (int index = 0; index< size; index++)
    result = result + V1[index]*V1[index];
  
  return result;
}

// Image with the sum square of the patch
void SS_Image(float *ima_in, float *ssi, int *neighborhoodsize, int *vol_size)
{
 
  int size = (2*neighborhoodsize[0]+1)*(2*neighborhoodsize[1]+1)*(2*neighborhoodsize[2]+1);
  float V1[size];
  
  for (int k=0; k < vol_size[2] ; k++)
  {
    for (int j=0; j < vol_size[1] ; j++)
    {
      for (int i=0; i < vol_size[0] ; i++)
      {
        Neiborghood(ima_in,i,j,k, neighborhoodsize, V1, vol_size, -1); // -1 to avoid normalization
        ssi[k*(vol_size[0]*vol_size[1])+(j*vol_size[0])+i] = L2_norm(V1, neighborhoodsize);
        //set_volume_real_value(ssi,i,j,k,0,0, L2_norm(V1, neighborhoodsize));
      }
    }
  }

}


// Patch extraction for FFT convolution
 void Patch_FFT(Volume ima_in, Volume patch, int x,int y,int z,int *neighborhoodsize, int *vol_size)
{
  int x_pos,y_pos,z_pos;
  bool is_outside;

  int nb_inside = 0;
  int count = 0;

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
        
         //Inversion of patch for correlation in transform domain
        
        if (is_outside) set_volume_real_value(patch, 2*neighborhoodsize[0]-a, 2*neighborhoodsize[1] - b, 2*neighborhoodsize[2] - c,0,0,0.0);
        
        else  set_volume_real_value(patch, 2*neighborhoodsize[0] - a, 2*neighborhoodsize[1] - b, 2*neighborhoodsize[2] - c,0,0, get_volume_real_value(ima_in, x_pos, y_pos, z_pos, 0, 0));

        
      }
    }
  }

}

// Search area extraction for FFT convolution
 void SearchArea_FFT(Volume ima_in, Volume search_area, int x,int y,int z,int *searchsize,int *neighborhoodsize, int *vol_size)
{
  int x_pos,y_pos,z_pos;
  bool is_outside;

  int nb_inside = 0;
  int count = 0;

  for (int c = 0; c<(2*searchsize[2] + 2*neighborhoodsize[2] +1);c++){
    for (int b = 0; b<(2*searchsize[1] + 2*neighborhoodsize[1] +1);b++){
      for (int a = 0;a<(2*searchsize[0] + 2*neighborhoodsize[0] +1);a++){

        is_outside = false;
        x_pos = x + a - searchsize[0] - neighborhoodsize[0];
        y_pos = y + b - searchsize[1] - neighborhoodsize[1];
        z_pos = z + c - searchsize[2] - neighborhoodsize[2];

        if ((z_pos < 0) || (z_pos > vol_size[2]-1)) is_outside = true;
        if ((y_pos < 0) || (y_pos > vol_size[1]-1)) is_outside = true;
        if ((x_pos < 0) || (x_pos > vol_size[0]-1)) is_outside = true;
        
        if (is_outside) set_volume_real_value(search_area,a,b,c,0,0,0.0);
        else set_volume_real_value(search_area,a,b,c,0,0,get_volume_real_value(ima_in, x_pos, y_pos, z_pos, 0, 0));
        
      }
    }
  }

}


// Convolution in transform domain
 void Convolution_FFT(Volume searcharea, Volume patch, Volume FFT_multi, int *FFTsize, pthread_mutex_t *mp)
{

  Volume tmp_p;
  Volume tmp_sa;

  float Re = 0.0;
  float Im = 0.0;
  
  prep_volume(&patch, &tmp_p);
  prep_volume(&searcharea, &tmp_sa);
  
  fft_volume(tmp_p, FALSE, 3, 0, mp);
  fft_volume(tmp_sa, FALSE, 3, 0, mp);
  
  for (int c = 0; c<(FFTsize[2]);c++)
  {
    for (int b = 0; b<(FFTsize[1]);b++)
    {
      for (int a = 0;a<(FFTsize[0]);a++)
      {
        
        Re = get_volume_real_value(tmp_sa, a, b, c, 0, 0) * get_volume_real_value(tmp_p, a, b, c, 0, 0) - get_volume_real_value(tmp_sa, a, b, c, 1, 0) * get_volume_real_value(tmp_p, a, b, c, 1, 0);
        
        Im = get_volume_real_value(tmp_sa, a, b, c, 1, 0) * get_volume_real_value(tmp_p, a, b, c, 0, 0) + get_volume_real_value(tmp_sa, a, b, c, 0, 0) * get_volume_real_value(tmp_p, a, b, c, 1, 0);
        
        
        set_volume_real_value(FFT_multi,a,b,c, 0 , 0, Re);
        set_volume_real_value(FFT_multi,a,b,c, 1 , 0, Im);
         
      }
    }
  }
  
  fft_volume(FFT_multi, TRUE, 3, 0, mp);

  delete_volume(tmp_p);
  delete_volume(tmp_sa);

}


//Main function
 void *Sub_denoise_fft_mt(void *arguments)
{

  nl_mean_fft_mt arg;
  arg=*(nl_mean_fft_mt *)arguments;
  int PID=0;
  PID=arg.thread_num;
  int ret;
  
  
  /* Creation of the required volumes */
  
  Volume patch,searcharea;
  char    *axis_order_fft[3] = {MIxspace, MIyspace, MIzspace};
  Real start[3] = {0,0,0};
  Real steps[3] = {1,1,1};
  
  patch = create_volume(3, axis_order_fft, NC_FLOAT, FALSE, 0.0, 0.0);
  set_volume_sizes(patch, arg.FFTsize);
  set_volume_starts(patch, start);
  set_volume_separations(patch, steps);
  alloc_volume_data(patch);
  for (int k = 0; k < arg.FFTsize[2]; k++)
    for (int j = 0; j < arg.FFTsize[1]; j++)    
        for (int i = 0; i < arg.FFTsize[0]; i++)
        set_volume_real_value(patch, i, j, k, 0, 0, 0.0);
  
  searcharea = create_volume(3, axis_order_fft, NC_FLOAT, FALSE, 0.0, 0.0);
  set_volume_sizes(searcharea, arg.FFTsize);
  set_volume_starts(searcharea, start);
  set_volume_separations(searcharea, steps);
  alloc_volume_data(searcharea);
  for (int k = 0; k < arg.FFTsize[2]; k++)
    for (int j = 0; j < arg.FFTsize[1]; j++)
        for (int i = 0; i < arg.FFTsize[0]; i++)
        set_volume_real_value(searcharea, i, j, k, 0, 0, 0.0);
  
  Volume FFT_multi;
  prep_volume(&searcharea, &FFT_multi); //storage of result
    
  float dist =0.0;
  double epsilon= 0.0001;

  for (int k=arg.debut; k < arg.fin ; k++)
  {
    double count = 0;
   
    for (int j=0; j < arg.vol_size[1] ; j++)
    {
     
      for (int i=0; i < arg.vol_size[0] ; i++)
      {
        float global_sum = 0.0;
        float average    = 0.0;
        float weight     = 0.0;
        float dist_min   = 1000000.0; 
                  
        Patch_FFT(arg.in,patch, i, j, k, arg.neighborhoodsize, arg.vol_size);
        SearchArea_FFT(arg.in,searcharea, i, j, k, arg.search, arg.neighborhoodsize, arg.vol_size);
        Convolution_FFT(searcharea, patch, FFT_multi, arg.FFTsize, arg.mp);
        
        for (int z = arg.neighborhoodsize[2]; z < arg.FFTsize[2] - (2*arg.neighborhoodsize[2] +1) ; z++)
        {
          for (int y = arg.neighborhoodsize[1]; y < arg.FFTsize[1] - (2*arg.neighborhoodsize[1] +1) ; y++)
          {
            for (int x = arg.neighborhoodsize[0]; x < arg.FFTsize[0] -  (2*arg.neighborhoodsize[0] +1) ; x++)
            {
        
              if ( 
                   (i + x - arg.neighborhoodsize[0] - arg.search[0]) >= 0 &&
                   (j + y - arg.neighborhoodsize[1] - arg.search[1]) >= 0 &&
                   (k + z - arg.neighborhoodsize[2] - arg.search[2]) >= 0 &&
                   (i + x - arg.neighborhoodsize[0] - arg.search[0]) < arg.vol_size[0] &&
                   (j + y - arg.neighborhoodsize[1] - arg.search[1]) < arg.vol_size[1] &&
                   (k + z - arg.neighborhoodsize[2] - arg.search[2]) < arg.vol_size[2] &&
                    x != arg.neighborhoodsize[0] + arg.search[0] &&
                    y != arg.neighborhoodsize[1] + arg.search[1] &&
                    z != arg.neighborhoodsize[2] + arg.search[2] 
                 )
                
              {
                
             
              dist = get_volume_real_value(arg.ssi, i, j, k, 0, 0)  +  
                      get_volume_real_value(arg.ssi, i + x - arg.neighborhoodsize[0] - arg.search[0], j + y - arg.neighborhoodsize[1] - arg.search[1], k + z - arg.neighborhoodsize[2]- arg.search[2] , 0, 0)
                  - 2 * get_volume_real_value(FFT_multi, x + (arg.neighborhoodsize[0] +1), y + (arg.neighborhoodsize[1] +1), z + (arg.neighborhoodsize[2] +1), 0, 0);

              if (dist < dist_min) dist_min = dist;
              weight = exp(- dist / (arg.filtering_param * arg.filtering_param));
              //weight = (1 / dist);
              global_sum = global_sum + weight;
              average = average + weight * get_volume_real_value(searcharea, x,  y,  z, 0, 0);

              }
              
           
            }
          }
        }
        
        /* Max weight is asigned to the central voxel */ 
        weight = exp(- dist_min / (arg.filtering_param * arg.filtering_param));
        global_sum = global_sum + weight;
        average = average + weight * get_volume_real_value(searcharea, arg.neighborhoodsize[0] + arg.search[0],  arg.neighborhoodsize[1] + arg.search[1], arg.neighborhoodsize[2] + arg.search[2] , 0, 0);
        
        if (global_sum != 0.0)
          set_volume_real_value(arg.out,i,j,k,0,0, (average/global_sum));
        else
          set_volume_real_value(arg.out,i,j,k,0,0,get_volume_real_value(arg.in,i,j,k,0,0));  
        
        
        }
      }
          

    if(debug)
      std::cout<<"Thread ( "<< std::setw(2) <<PID+1<< " ) has finished the slice :"<<std::setw(3)<< k <<std::endl;

//     ret = pthread_mutex_lock(arg.mp);
// //    *arg.mean_neighboors = *arg.mean_neighboors + count/(arg.vol_size[0]*arg.vol_size[1]);
//     ret = pthread_mutex_unlock(arg.mp);

 
  }
  
  delete_volume(patch);
  delete_volume(searcharea);

  if(debug)
    std::cout<<"End of the thread : "<<std::setw(2)<<PID+1<<std::endl;
  pthread_exit(0);


}



//Creation of the threads
 void denoise_mt_fft( Volume in,Volume out,Volume mean_map, Volume var_map, double filtering_param, int *neighborhoodsize, int *search,int testmean,int testvar,double m_min,double v_min,int weight_method, int * vol_size)
{

  int ii;
  pthread_t	*thread;
  nl_mean_fft_mt	*arguments;
  void	*retval;
  thread = (pthread_t *) calloc(nb_thread, sizeof(pthread_t));
  arguments=(nl_mean_fft_mt *)calloc(nb_thread,sizeof(nl_mean_fft_mt));

  float mean_neighboors = 0;
  
  int powsize = 1;
  
  // We search the 2^n size of the volume to convoluate
  while ( (((2*search[0] + 1) + (2*neighborhoodsize[0]+1) ) > pow(2,powsize)) || (((2*search[1] + 1) + (2*neighborhoodsize[1]+1)) > pow(2,powsize)) || (((2*search[2] + 1) + (2*neighborhoodsize[2]+1)) > pow(2,powsize)) )
    powsize ++;
         
  //std::cout<<"powsize ="<<pow(2,powsize)<<std::endl;
  int FFTsize[3] = {(int)pow(2,powsize) , (int)pow(2,powsize) , (int)pow(2,powsize)};

  int ret;
  pthread_mutex_t mp = PTHREAD_MUTEX_INITIALIZER;
  /* initialize a mutex to its default value */
  ret = pthread_mutex_init(&mp, NULL);
 
  
  Volume ssi;
  ssi = copy_volume_definition(in, NC_FLOAT, FALSE, 0.0, 0.0);
  for (int i = 0; i < vol_size[0]; i++)
    for (int j = 0; j < vol_size[1]; j++)
      for (int k = 0; k < vol_size[2]; k++)
        set_volume_real_value(ssi, i, j, k, 0, 0, 0.0);
  
  SS_Image(in, ssi, neighborhoodsize, vol_size); 


  for(ii=0;ii<nb_thread;ii++)
  {
    arguments[ii].in=in;
    arguments[ii].out=out;
    arguments[ii].filtering_param=filtering_param;
    arguments[ii].search=search;
    arguments[ii].vol_size=vol_size;
    arguments[ii].FFTsize=FFTsize;
    arguments[ii].neighborhoodsize=neighborhoodsize;
    arguments[ii].ssi=ssi;
    arguments[ii].thread_num=ii;
    arguments[ii].mp=&mp;
    arguments[ii].debut=(int)(ii*vol_size[2]/nb_thread);
    arguments[ii].fin=(int)((ii+1)*vol_size[2]/nb_thread);
  }


  if(debug) {
    std::cout<<"\n Creation of threads\n"<<std::endl;
    std::cout<<" Number of slices : "<<vol_size[2]<<"\n"<<std::endl;
  }
  
  for(ii=0;ii<nb_thread;ii++)
    if(pthread_create(&thread[ii],NULL,Sub_denoise_fft_mt,&arguments[ii]))
    {
      std::cerr<<"Creation de thread impossible"<<std::endl;
      exit(1);
    }

  for(ii=0;ii<nb_thread;ii++)
    if(pthread_join(thread[ii],&retval))
    {
      std::cerr<<"Synchronisation de thread impossible"<<std::endl;
      exit(1);
    }
    
    
  free(thread);
  free(arguments);

  if(verbose) {
    std::cout << "\nMean number of neighboors per voxel: " << mean_neighboors/vol_size[2] << std::endl;
    
    std::cout << "\nVolume denoising is finished"<< std::endl;
  }
}
