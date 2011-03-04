#ifndef NL_MEANS_BLOCK_HH
#define NL_MEANS_BLOCK_HH

#include "nl_means.h"
#include "nl_means_utils.h"

using namespace std;

// Structur used for multithreading programming
struct nl_mean_block_mt {
  volatile float *Estimate;
  volatile float *Label;
  float *in;
  float *mean_map;
  float *var_map;
  double filtering_param;
  double beta;
  int debut;
  int fin;
  int *neighborhoodsize;
  int *search;
  int *vol_size;
  int testmean;
  int testvar;
  double m_min;
  double v_min;
  int weight_method;
  int thread_num;
  int b_space;
  pthread_mutex_t *mp;
  float *mean_neighboors;
  float *hallucinate;
};

// Function which compute the weighted average for one block
 void Average_block(float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float weight, int* vol_size);

// Function which compute the weighted average for one block
 void Average_block_Rician(float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float weight, int* vol_size);


// Function which computes the value assigned to each voxel
 void Value_block(pthread_mutex_t *mp, volatile float *Estimate, volatile float* Label, float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float global_sum, int* vol_size);

// Function which computes the value assigned to each voxel
 void Value_block_Rician(pthread_mutex_t *mp, volatile float *Estimate, volatile float* Label, float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float global_sum, int* vol_size, float filtering_param);

//Main function
void *Sub_denoise_block_mt(void *arguments);

//Creation of the threads
void denoise_block_mt(float *in, float *out, float *mean_map, float *var_map, double filtering_param, double beta, int *neighborhoodsize, int *search,int testmean,int testvar,double m_min, double v_min,int weight_method, int b_space, int *vol_size, float *hallucinate);

#endif
