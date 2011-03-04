#ifndef NL_MEANS_H
#define NL_MEANS_H

#include <iostream> 
#include <iomanip>
#include <pthread.h>
#include <volume_io.h>

extern int testmean,testvar, nb_thread, block;

// Structure used for multithreading programming
struct nl_mean_mt {
  float  *in;
  volatile float *out;
  float *mean_map;
  float *var_map;
  float *weight;
  float *w_max;
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
  pthread_mutex_t *mp;
  float *mean_neighboors;
  float *hallucinate;
};

//Computation of the L2 norm between 2 neighborhoods
 float L2_norm(float *V1, float *V2, int *neighborhoodsize);

//Computation of the Pearson ditance between two neighborhoods (for US image)
 float Pearson_distance(float *V1, float *V2, int *neighborhoodsize);

//Weighting function
 float Weight(float dist,float beta,float h);

// Neighborhoods extraction
 void Neiborghood(float *ima_in,int x,int y,int z,int *neighborhoodsize,float *neighborhoodvec, int *vol_size, int weight_method);

//Main function
 void *Sub_denoise_mt(void *arguments);

//Creation of the threads
 void denoise_mt(float *in,float *out, float *mean_map, float *var_map, double filtering_param,double beta, int *neighborhoodsize, int *search,int testmean,int testvar,double m_min,double v_min,int weight_method, int * vol_size,float *hallucinate);

#endif
