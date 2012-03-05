// Parzen window density estimation routines. 
// MRFSEG + GAMIXTURE software bundle: Voxel classification in volumetric 
// images using genetic algorithms and Markov random fields.
// The genetic algorithm is described in: 
// [1]  J. Tohka , E. Krestyannikov, I.D. Dinov , A. MacKenzie-Graham,
// D.W. Shattuck , U. Ruotsalainen, and A.W. Toga. 
// Genetic algorithms for finite mixture model based voxel classification
// in neuroimaging.  
// IEEE Transactions on Medical Imaging  , 26(5):696 - 711, 2007.
// (C) Jussi Tohka, 2005 - 2007 
// Laboratory of NeuroImaging, Department of Neurology, UCLA Medical School, CA, USA
// Institute of Signal Processing, Tampere University of Technology, Finland  
// Version 1.1
// e-mail jussi.tohka@tut.fi
// *****************************************************
// ******************************************************
// Permission to use, copy, modify, and distribute this software 
// for any purpose and without fee is hereby
// granted, provided that the above copyright notice appears in all
// copies.  The author, Tampere University of Technology and University of
// California, Los Angeles make no representations
// about the suitability of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
// *****************************************************
// Version 1.1 offers some speed ups compared to the version 1.0

#include "parzen.h"
#include <iostream>
#include <fstream>

//todo: use std instead
#ifdef   MAX
#undef   MAX
#endif
#define  MAX( x, y )  ( ((x) >= (y)) ? (x) : (y) )

#ifdef   MIN
#undef   MIN
#endif
#define  MIN( x, y )  ( ((x) <= (y)) ? (x) : (y) )


#define ROUND(x) (ceil((x) - 0.5))
#define ARE_EQUALF(x,y) (fabs((x) - (y)) < 0.0000001)   
#define FLOAT_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

using namespace std;

// computes the x coordinates for the estimate
void computeX(PdfEstimate* hatf,const FloatImage& img, const LabelImage& mask,float percentage )
{
  float maximum,minimum,step;
  int i;  

  if(fabs(percentage - 1) < 0.00001) { 
    maximum = img.imageMax(mask);
  }
  else {
    maximum = img.imageKth(mask,percentage);
  }
  minimum = img.imageMin(mask);
  hatf->x.resize(hatf->n);
  step = (maximum - minimum)/((float) (hatf->n - 1));
  for(i = 0;i < hatf->n;i++) {
    hatf->x[i] = minimum +(i * step);
  }
}

void computeY(PdfEstimate* hatf,const FloatImage& img,const FloatImage& mask)
{
  int x,y,z,imgsize,j,enditer,startiter;
  float masksum;  
  float c;
  float beta;
  float step,ysum;

  hatf->y.resize(hatf->n,0);
 
  for(j = 0;j < hatf->n;j++) {
    hatf->y[j] = 0;
  }   

  step = hatf->x[1] - hatf->x[0]; // assumes equal x-axis spacing 
  beta = 2 * (hatf->sigma) * (hatf->sigma);
  c = sqrt(2 * M_PI) * (hatf->sigma);
  masksum = 0;
  ysum = 0;  

  // imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim);
  for(x = 0;x < img.header.x_dim;x++) {
    for(y = 0;y < img.header.y_dim;y++) {
      for(z = 0;z < img.header.z_dim;z++) {
       
        if(mask.getVoxelValue(x,y,z) > 0.001) {
          masksum = masksum + mask.getVoxelValue(x,y,z);
          startiter = (int) floor((img.getVoxelValue(x,y,z) - 5*(hatf->sigma) - hatf->x[0])/step); // speeding up. Assumes equal
          enditer = (int) ceil((img.getVoxelValue(x,y,z) + 5*(hatf->sigma) - hatf->x[0])/step) + 1;   // x-axis spacing 
          startiter = MAX(startiter,0);
          enditer = MIN(enditer,hatf->n); 
          for(j = startiter;j < enditer;j++) {
            hatf->y[j] = hatf->y[j] +mask.getVoxelValue(x,y,z)*exp(-(pow((hatf->x[j] - img.getVoxelValue(x,y,z)),2)/beta));
          }
        }
      }
    }
  }
  
  for(j = 0;j < hatf->n;j++) {
    hatf->y[j] = hatf->y[j]/(c*masksum);
    ysum = ysum + hatf->y[j];
  } 
  ysum = ysum*step;
//  cout << "Parzen integrand " << ysum << endl;  
  // take care that the pdf intergrates to 1. Note that we study intergration range from x[0] - step/2 to x[n] + step/2 
  // This precaution is because the speed up technique
  for(j = 0;j < hatf->n;j++) {
    hatf->y[j] = hatf->y[j]/(ysum);
  } 
 
}

// writes pdf estimate to a text file
// format of the text file:
// 1st line : float sigma int n 
// 2nd - nth + 1 line: x[i] y[i] 

bool writeEstimate(char* filename,PdfEstimate* hatf, bool overwrite) 
{
  int i;
  
  if(! overwrite) {
    ifstream ifile(filename);
    if(ifile) {
      ifile.close();
      return(false);
    }
  }
  ofstream ofile(filename);
  if(!ofile) return(false);
  ofile << hatf->sigma << " " << hatf->n << "\n";
  for(i = 0;i < hatf->n;i++) {
    ofile << hatf->x[i] << " " << hatf->y[i] << "\n"; 
  }
  ofile.close();
  return(true);
}

bool readEstimate(char* filename,PdfEstimate* hatf)
{
  int i;

  ifstream ifile(filename);
  if(!ifile) return(false);
  ifile >> hatf->sigma;
  ifile >> hatf->n;
  hatf->x.resize(hatf->n);
  hatf->y.resize(hatf->n);
  for(i = 0;i < hatf->n;i++) {
    ifile >> hatf->x[i];
    ifile >> hatf->y[i];
  }
  return(true);
}
