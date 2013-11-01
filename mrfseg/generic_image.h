
// ******************************************************
// Functions to manipulate Analyze 7.5 images 
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
// Thanks to David Shattuck with help with these functions. 
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

// (C) Jussi Tohka, 2005  (Based partly on the code provided by David Shattuck)
// Laboratory of NeuroImaging, Department of Neurology, UCLA Medical School, CA, USA
// Institute of Signal Processing, Tampere University of Technology, Finland  
// Version 1.0, for internal distribution only.
// e-mail jussi.tohka@{tut.fi,loni.ucla.edu}
// *****************************************************

#ifndef __GENERIC_IMAGE_H__
#define __GENERIC_IMAGE_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H


#include <vector>
#include <algorithm>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_VALUES_H
#include <values.h>
#else
#include <limits.h>
#include <float.h>
#define MAXSHORT   SHRT_MAX
#define MAXINT     INT_MAX 
#define MAXDOUBLE  DBL_MAX
#define MAXFLOAT   FLT_MAX

#define MINSHORT   SHRT_MIN
#define MININT     INT_MIN 
#define MINDOUBLE  DBL_MIN
#define MINFLOAT   FLT_MIN

#endif 

#include <math.h>


struct ImageHeader {
	ImageHeader():start(3,0),step(3,0),dims(3),dir_cos(3,std::vector<double>(3,0.0))
	{
    //Identity
    dir_cos[0][0]=dir_cos[1][1]=dir_cos[2][2]=1.0;
	}
	int	dims;		  
	int	x_dim;		
	int	y_dim;		
	int	z_dim;		
  std::vector<double> start;
  std::vector<double> step;
  std::vector<std::vector<double> > dir_cos;
};


class FloatImage;

//todo: this maybe trivially parametrized as template
class LabelImage {
  public:
  ImageHeader header;
  std::vector<unsigned char> data; // If memory is valuable
  
  unsigned char getLabelValue(int x, int y, int z) const
  {
    return data[x + y*header.x_dim + z*header.y_dim*header.x_dim];
  };

  // this is the same as the previous one exept returns 0 if we are ouside the image

  unsigned char getSafeLabelValue(int x, int y, int z) const
  {
    if((x < 0) || (y < 0) || (z < 0) || (x > (header.x_dim - 1)) || 
       (y > (header.y_dim - 1)) || (z > (header.z_dim - 1))) 
      return(0);
    else  return data[x + y*(header.x_dim) + z*(header.y_dim)*(header.x_dim)];
  }

  void putLabelValue(int x, int y, int z, unsigned char val) 
  {
    data[x + y*(header.x_dim) + z*(header.y_dim)*(header.x_dim)] = val;
  };
  
  unsigned char labelMax(void);
  bool copyImage(const LabelImage& source);
  bool copyImage(const FloatImage& source);
  bool newImage(const LabelImage& alike);  
  bool newImage(const FloatImage& alike);
  void thresholdImage(const FloatImage& source,float val);
  void erode3D(void);
};

//todo: this maybe trivially paramtrized as template
class FloatImage {
  public:
  ImageHeader header;
  std::vector<float> data; // Internal data presented using floats 
  
  float getVoxelValue(int x, int y, int z) const
  {
    return data[x + y*(header.x_dim) + z*(header.y_dim)*(header.x_dim)];
  };
  
  float getSafeVoxelValue(int x, int y, int z) const
  {
    if((x < 0) || (y < 0) || (z < 0) || (x > (header.x_dim - 1)) || 
       (y > (header.y_dim - 1)) || (z > (header.z_dim - 1))) 
      return(0);
    else  return data[x + y*(header.x_dim) + z*(header.y_dim)*(header.x_dim)];
  }

  void putVoxelValue(int x, int y, int z, float val) 
  {
    data[x + y*(header.x_dim) + z*(header.y_dim)*(header.x_dim)] = val;
  };

  void thresholdImage(float val);
  
  float imageMax(const LabelImage& mask) const;
  float imageMin(const LabelImage& mask) const;
  
  bool copyImage(const FloatImage& source);
  bool copyImage(const LabelImage& source);
  bool newImage(const FloatImage& alike);  
  bool findClosestNonZero(int x,int y, int z,int* cx, int* cy, int* cz,float tr) const;
  float imageKth(const LabelImage& mask, float percentage) const ;

  void medianFilter3D(const LabelImage& mask, int windowX, int windowY, int windowZ);
  void knnFilter(const LabelImage& mask, int windowLenX, int windowLenY, int windowLenZ, int knn);
}; 

bool binaryMask(const FloatImage& img,const LabelImage& mask,float tr); 
void averageFilterZeros(FloatImage& img,const FloatImage& refImg,float tr,int times);

  
inline char maxArg(float* pval,int n)
{
  float maximum;
  int i,index;
  
  maximum = pval[0];
  index = 0;
  for(i = 1;i < n;i++) {
    if(pval[i] > maximum) {
      index = i;
      maximum = pval[i];
    }
  }
  return( (char) index);
}
#endif //GENERIC_IMAGE

