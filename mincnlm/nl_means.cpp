#include "nl_means.h"

extern int      verbose;
extern int      debug;

//Computation of the L2 norm between 2 neighborhoods
 float L2_norm(float *V1, float *V2, int *neighborhoodsize)
{

  double result,difference;
  result = 0;
  int size = (2*neighborhoodsize[0]+1)*(2*neighborhoodsize[1]+1)*(2*neighborhoodsize[2]+1);
  for (int index = 0; index< size; index++) {
    difference = (V1[index]-V2[index]);
    result = result + difference*difference;
  }

  return result;
}

//Computation of the Pearson ditance between two neighborhoods (for US image)
 float Pearson_distance(float *V1, float *V2, int *neighborhoodsize)
{

  float result,difference;
  result = 0;
  int size = (2*neighborhoodsize[0]+1)*(2*neighborhoodsize[1]+1)*(2*neighborhoodsize[2]+1);
  for (int index = 0; index< size; index++) {
    difference = (float) (V1[index]-V2[index]);

   //avoid division by 0 + symetric distance
   result = result + ((difference*difference )/ ((V1[index] + V2[index])/2 + 0.0001 )); 
  }

  return result;
}


//Weighting function
 float Weight(float dist,float beta,float h)
{
  if (h != 0)
    return exp((-dist/(2.0*beta*h*h)));
  else
    return 0;
}


// Neighborhoods extraction
 void Neiborghood(float *ima_in,int x,int y,int z,int *neighborhoodsize,float *neighborhoodvec, int *vol_size, int weight_method)
{
  int x_pos,y_pos,z_pos;
  bool is_outside;

  int nb_inside = 0;
  int count = 0;

  int size = (2*neighborhoodsize[0]+1)*(2*neighborhoodsize[1]+1)*(2*neighborhoodsize[2]+1);

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
        if (is_outside) neighborhoodvec[count] = 0;
        else {
          
          neighborhoodvec[count] =  ima_in[z_pos*(vol_size[0]*vol_size[1])+(y_pos*vol_size[0])+ x_pos];
          nb_inside++;
        }
        count++;
      }
    }
  }
  
  if ((weight_method == 0) || (weight_method == 2))
    {
    
      for (int index=0; index <size;index++)
      {
        // Normalization by the number of elements 
        neighborhoodvec[index] = (neighborhoodvec[index]/sqrt((float) nb_inside));
      }
    }
}

//Main function
 void *Sub_denoise_mt(void *arguments)
{

  nl_mean_mt arg;
  arg=*(nl_mean_mt *)arguments;
  int PID=0;
  PID=arg.thread_num;
  int ret;
  int x_min,x_max,y_min,y_max, z_min, z_max;
  int NbElement2 = (2*arg.neighborhoodsize[0]+1)*(2*arg.neighborhoodsize[1]+1)*(2*arg.neighborhoodsize[2]+1);
  
  typedef float Voisinages[NbElement2];
  Voisinages V1,V2;
  
  float w_max = 0;

  double epsilon= 0.0001;
  
  int offset_k = 0;
  int offset_j = 0;
  int offset_kk = 0;
  int offset_jj = 0;

  for (int k=arg.debut; k < arg.fin ; k++)
  {
    double count = 0;
    offset_k = k*(arg.vol_size[0]*arg.vol_size[1]);
   
    for (int j=0; j < arg.vol_size[1] ; j++)
    {
     
      offset_j = j*arg.vol_size[0];
      
      for (int i=0; i < arg.vol_size[0] ; i++)
      {
        float global_sum = 0;
        float average    = 0;
      


        x_min = MAX(0,i-arg.search[0]);      x_max = MIN(arg.vol_size[0]-1,i+arg.search[0]);
        y_min = MAX(0,j-arg.search[1]);      y_max = MIN(arg.vol_size[1]-1,j+arg.search[1]);
        z_min = MAX(0,k-arg.search[2]);      z_max = MIN(arg.vol_size[2]-1,k+arg.search[2]);

         // Only the voxel with a local mean and variance different to zero are processed. 
        // By this way the computational time is reduced since the background is not taken into account (or the B-scan mask for US image)
       
       
        if ((arg.mean_map[offset_k + offset_j + i] > epsilon) && (arg.var_map[offset_k + offset_j + i] > epsilon))
        {

          if ((arg.weight_method == 0) || (arg.weight_method == 2)) 
            Neiborghood(arg.in,i,j,k,arg.neighborhoodsize,V1,arg.vol_size,arg.weight_method);
//             
          if (arg.weight_method == 1)
            Neiborghood(arg.in,i,j,k,arg.neighborhoodsize,V1,arg.vol_size,arg.weight_method);

          w_max = 0;
          
          for (int kk = z_min; kk <= z_max; kk++)
          {
            offset_kk = kk*(arg.vol_size[0]*arg.vol_size[1]);
            
            for (int jj = y_min; jj <= y_max; jj++)
            {
              
              offset_jj = jj*arg.vol_size[0];
              
              for (int ii = x_min; ii <= x_max; ii++)
              {
		
                // To avoid the compute the weigh with null patch and the division by zero in ratio and ratio 2
                
                if ((arg.mean_map[offset_kk + offset_jj + ii] > epsilon) && (arg.var_map[offset_kk + offset_jj + ii] > epsilon))   
                {

                  
                  float ratio  = arg.mean_map[offset_k + offset_j + i]/arg.mean_map[offset_kk + offset_jj + ii];

                  float ratio2 = arg.var_map[offset_k + offset_j + i]/arg.var_map[offset_kk + offset_jj + ii];
                  
                  float weight = 0;

			
                      if ((arg.weight_method == 0) || (arg.weight_method == 2))
                      {
                            
                          if (   (arg.m_min <= ratio) && (ratio <= (1/arg.m_min)) 
                              && (arg.v_min <= ratio2) && (ratio2 <= (1/arg.v_min)) )
                          {
                            
                            
                            if ((ii != i) || (jj != j) || (kk != k))
                            {
                              float val_in;
                              count = count + 1;
                              Neiborghood(arg.in,ii,jj,kk,arg.neighborhoodsize,V2,arg.vol_size,arg.weight_method);
                              
                              weight = Weight(L2_norm(V1,V2,arg.neighborhoodsize),arg.beta,arg.filtering_param);
                              global_sum = global_sum + weight;
                              
                              
                              val_in=arg.hallucinate?arg.hallucinate[offset_kk + offset_jj + ii]:arg.in[offset_kk + offset_jj + ii];//VF
                              
                              //Rician denoising
                              if (arg.weight_method == 2) 
                                average = average + val_in * val_in *weight;
                              else
                                average = average + val_in * weight;
                              
                              if (weight > w_max) w_max = weight;
                            }
                          }
                        }
                        
                       else
                        {
                            
                         if (  (arg.m_min <= ratio) && (ratio <= (1/arg.m_min)))
                          {
                            if ((ii != i) || (jj != j) || (kk != k))
                            {
                              float val_in;
                              count = count + 1;
                              //Pearson distance computation
                              weight = Weight(Pearson_distance(V1,V2,arg.neighborhoodsize),arg.beta,arg.filtering_param);
                              global_sum = global_sum + weight;
                              val_in=arg.hallucinate?arg.hallucinate[offset_kk + offset_jj + ii]:arg.in[offset_kk + offset_jj + ii]; //VF
                              average = average + val_in * weight;
                              if (weight > w_max) w_max = weight;

                             }
                           }
                         }
                }
              }
            }
          }
	
       
          
          global_sum = global_sum + w_max;
          
          double denoised_value = 0.0;
          float val_in=arg.hallucinate?arg.hallucinate[offset_k + offset_j + i]:arg.in[offset_k + offset_j + i];//VF
         
        
          
          
          // Rician denoising
          if (arg.weight_method == 2)
          {
            average = average + val_in * val_in * w_max;
            
            if (global_sum != 0.0) 
            {  
             
              denoised_value = (average/global_sum) - (2*arg.filtering_param * arg.filtering_param); 
              
              if (denoised_value > 0.0)
                denoised_value = sqrt (denoised_value);
              else 
                denoised_value = 0.0;
              
              arg.out[offset_k + offset_j + i] = denoised_value;
            }
            else
            {
              arg.out[offset_k + offset_j + i] =  val_in;
            }
          }
          // Gaussian and Speckle denoising
          else
          {
           
            average = average + val_in*w_max; 
          
            if (global_sum != 0.0)
              arg.out[offset_k + offset_j + i] = average/global_sum;
            else
              arg.out[offset_k + offset_j + i] =  val_in;
          }
           
	
        }
      }
    }

    if(debug)
      std::cout<<"Thread ( "<< std::setw(2) <<PID+1<< " ) has finished the slice :"<<std::setw(3)<< k <<std::endl;

    ret = pthread_mutex_lock(arg.mp);
    *arg.mean_neighboors = *arg.mean_neighboors + count/(arg.vol_size[0]*arg.vol_size[1]);
    ret = pthread_mutex_unlock(arg.mp);
  }
  if(debug)
    std::cout<<"End of the thread : "<<std::setw(2)<<PID+1<<std::endl;
  
  pthread_exit(0);
}


//Creation of the threads
 void denoise_mt(float *in, float *out, float *mean_map, float *var_map, double filtering_param, double beta, int *neighborhoodsize, int *search,int testmean,int testvar,double m_min,double v_min,int weight_method, int * vol_size, float *hallucinate)
{

  int ii;
  pthread_t	*thread;
  nl_mean_mt	*arguments;
  void	*retval;
  thread = (pthread_t *) calloc(nb_thread, sizeof(pthread_t));
  arguments=(nl_mean_mt *)calloc(nb_thread,sizeof(nl_mean_mt));

  float mean_neighboors = 0;

  int ret;
  pthread_mutex_t mp = PTHREAD_MUTEX_INITIALIZER;
  /* initialize a mutex to its default value */
  ret = pthread_mutex_init(&mp, NULL);
  
  /* local and local variance storage */
  float *w_max, *weight, *storage;
  w_max = new float[vol_size[0]*vol_size[1]*vol_size[2]];
  weight = new float[vol_size[0]*vol_size[1]*vol_size[2]];

  for (int i = 0; i < vol_size[0] *  vol_size[1] * vol_size[2]; i++)
  {
    w_max[i] = 0.0;
    weight[i] = 0.0;
  
  }

  for(ii=0;ii<nb_thread;ii++)
  {
    arguments[ii].in=in;
    arguments[ii].out=out;
    arguments[ii].w_max=w_max;
    arguments[ii].weight=weight;
    arguments[ii].filtering_param=filtering_param;
    arguments[ii].beta=beta;
    arguments[ii].search=search;
    arguments[ii].vol_size=vol_size;
    arguments[ii].neighborhoodsize=neighborhoodsize;
    arguments[ii].mean_map=mean_map;
    arguments[ii].var_map=var_map;
    arguments[ii].testmean=testmean;
    arguments[ii].testvar=testvar;
    arguments[ii].m_min=m_min;
    arguments[ii].v_min=v_min;
    arguments[ii].weight_method = weight_method;
    arguments[ii].thread_num=ii;
    arguments[ii].mp=&mp;
    arguments[ii].mean_neighboors=&mean_neighboors;
    arguments[ii].debut=(int)(ii*vol_size[2]/nb_thread);
    arguments[ii].fin=(int)((ii+1)*vol_size[2]/nb_thread);
    arguments[ii].hallucinate=hallucinate;
  }

  if(debug) {
    std::cout<<"\n Creation of threads\n"<<std::endl;
    std::cout<<" Number of slices : "<<vol_size[2]<<"\n"<<std::endl;
  }
  for(ii=0;ii<nb_thread;ii++)
    if(pthread_create(&thread[ii],NULL,Sub_denoise_mt,&arguments[ii]))
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
