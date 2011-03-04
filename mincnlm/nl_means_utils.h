#ifndef NL_MEANS_UTILS_H
#define NL_MEANS_UTILS_H

#include "nl_means.h"

float FindMinVolume(float *ima_in, int *vol_size);
float FindMaxVolume(float *ima_in, int *vol_size);

// Linear scaling between 0 - 255
void LinearScaling(float *ima_in, double mini, double maxi, int *vol_size);

// Linear rescaling between mini - maxi
void LinearReScaling(float *ima_in, double mini, double maxi, int *vol_size);

float Neighborhood_Mean(float *ima_in,int x,int y,int z,int *neighborhoodsize, int * vol_size);

void Preprocessing(float *ima_in, float *mean_map,int *neighborhoodsize, int *vol_size);

float Neighborhood_Var(float  *ima_in, float *ima_mean, int x, int y, int z,int *neighborhoodsize, int * vol_size);

void Preprocessing2(float *ima_in, float *mean_map, float *var_map,int *neighborhoodsize, int * vol_size);

//STd estimation of the gaussian noise by pseudo residu
void Variance_Estimation_1(float *ima_in, int weight_method, double &s_constant, int * vol_size);

#endif
