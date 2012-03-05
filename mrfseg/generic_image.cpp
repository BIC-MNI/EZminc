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
// (C) Jussi Tohka, 2005 
// Laboratory of NeuroImaging, Department of Neurology, UCLA Medical School, CA, USA
// Institute of Signal Processing, Tampere University of Technology, Finland  
// Version 1.0, for internal distribution only.
// e-mail jussi.tohka@{tut.fi,loni.ucla.edu}
// *****************************************************

#include "generic_image.h"
#include <functional>
using namespace std;

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

bool FloatImage::copyImage(const FloatImage& source) 
{
  header = source.header;
  data = source.data;
  return true;
}

bool FloatImage::copyImage(const LabelImage& source) 
{
  header = source.header;
  data.resize(source.data.size());
  for(int i = 0;i < data.size();i++) {
    data[i] = source.data[i];
  }
  return true;
}

bool LabelImage::copyImage(const LabelImage& source)
{
  header = source.header;
  data = source.data;
  return true;
}

bool LabelImage::copyImage(const FloatImage& source)
{
  header = source.header;
  data.resize(source.data.size());
  for(int i = 0;i < data.size();i++) {
    data[i] = (unsigned char)source.data[i];
  }
  return true;
}

bool FloatImage::newImage(const FloatImage& alike) 
{
  header = alike.header;
  data.resize(alike.data.size(),0.0);
  return true;
}

bool LabelImage::newImage(const LabelImage& alike)
{
  header = alike.header;
  data.resize(alike.data.size(),0);
  return true;
}

bool LabelImage::newImage(const FloatImage& alike)
{
  header = alike.header;
  data.resize(alike.data.size(),0);
  return true;
}

void LabelImage::thresholdImage(const FloatImage& source,float tr)
{
  newImage(source);
  for(int i = 0;i < data.size();i++) {
    if(source.data[i] > tr) data[i] = 1;
    else data[i] = 0;
  }
}

// finds the closest voxel with the value greater than tr to the voxel (x,y,z)
// extremely slow but correct 
// returns false if no voxel value exceeds the threshold

bool FloatImage::findClosestNonZero(int x,int y, int z,int* cx, int* cy, int* cz,float tr) const
{
  int i,j,k;
  float d;
  bool found; 
  int windowSize = 11;
  // Initialize
  *cx = x;  // closest x coordinate
  *cy = y;  // closest y coordinate
  *cz = z;  // closest z coordinate
  found = false;
  // if value at (x,y,z) > tr then ok
  if(getVoxelValue(*cx,*cy,*cz) > tr) return(true);
  d = pow((float) (header.x_dim + header.y_dim  + header.z_dim),2);
  for(i = (-windowSize);i < (windowSize + 1);i++) {
    for(j = (-windowSize);j < (windowSize + 1) ;j++) {
      for(k = (-windowSize);k < (windowSize + 1);k++) {
         if(getSafeVoxelValue(x + i,y + j,z + k) > tr) {
           if( (pow((float) i,2) + pow((float) j,2) + pow((float) k,2)) < d) {
             d = (pow((float) i,2) + pow((float) j ,2) + pow((float) k,2));
             *cx = x + i;
             *cy = y + j;
             *cz = z + k;
             found = true;
           }
         }
       }
     }
   }
   return(found);
}

void averageFilterZeros(FloatImage& img,const FloatImage& refImg,float tr,int times)
{
  int x,y,z,i;
  int x1,y1,z1;
  double tmpf;  

  for(i = 0;i < times;i++) {
    for(x = 0;x < img.header.x_dim;x++) {
      for(y = 0;y < img.header.y_dim;y++) {
        for(z = 0;z < img.header.z_dim;z++) {
          if(refImg.getVoxelValue(x,y,z) < tr) {
            tmpf = 0.0;
            for(x1 = -3;x1 < 4;x1++) {
              for(y1 = -3;y1 < 4;y1++) {
                for(z1 = -3;z1 < 4;z1++) {
                  tmpf = tmpf + img.getSafeVoxelValue(x + x1,y + y1,z + z1);
                }
              }
            }
            img.putVoxelValue(x,y,z,tmpf/(7*7*7));
          }
        }
      }
    }
  } 
}

// finds a value such that exatcly a certain percentage of voxel values (in img) 
// are smaller than it. Considers only voxels within the mask.
// bubble sort inside the mask.... not very fast?
float FloatImage::imageKth(const LabelImage& mask, float percentage) const
{
  unsigned int i,j,k,l,m,n;
  
  n = 0;
  for(i = 0;i < data.size();i++) {
    if(mask.data[i] >=1 ) {
      n++;
    }
  }
  k = (unsigned int) floor(percentage*n);
  if(k > n) return(0);
  if(k < 0) return(0);
  std::vector<float> a(n);
  n = 0;
  for(i = 0;i < data.size();i++) {
    if(mask.data[i] >=1 ) {
      a[n] = data[i];
      n++;
    }
  }
 std::sort(a.begin(),a.end(),std::less<float>());  
 return a[k];
}

// Median filters image img; 
// Assumes windowLengths are odd integers
// for example (medianFilter3D(..., ..., 1, 1, 1);
// uses 3 x 3 x 3 mask.
//
//WARNING: this seem to be a recursive filter , is it what's expected!?
void FloatImage::medianFilter3D(const LabelImage& mask, int windowLenX, int windowLenY, int windowLenZ)
{
  int windowSize;
  int x,y,z,i,j,k,l,m;
  int n;
  int med;
  float f;

  windowSize = (2*windowLenX + 1)*(2*windowLenY + 1)*(2*windowLenZ + 1);  
  med = (windowSize - 1)/2;
  std::vector<float> window(windowSize);


  for(x = 0;x < header.x_dim;x++) {
    for(y = 0;y < header.y_dim;y++) {
      for(z = 0;z < header.z_dim;z++) {
        if(mask.getLabelValue(x,y,z) >= 1) {
          n = 0;
          for(i = (-windowLenX);i < (windowLenX + 1);i++) {
            for(j = (-windowLenY);j < (windowLenX + 1);j++) {
              for(k = (-windowLenZ);k < (windowLenX + 1);k++) {
                window[n] = getSafeVoxelValue(x + i, y + j, z + k);
                n++;
              }
            }
          }
          std::sort(window.begin(),window.end(),std::less<float>());
          putVoxelValue(x,y,z,window[med]);
        }
      }
    }
  }
}


void FloatImage::thresholdImage(float val) {
  for(int i = 0;i< data.size();i++) {
    if(data[i] > val) 
      data[i] = 1.0;
    else 
      data[i] = 0;
  }
}   

float FloatImage::imageMax(const LabelImage& mask)  const
{   
  int i,j,k;
  float maximum;
  
  maximum = MINFLOAT;
  //todo this is should be just straigth walk through 2 buffers, without nested loops!
  for(i = 0;i < header.x_dim;i++) {
    for(j = 0;j < header.y_dim;j++) {
       for(k = 0;k < header.z_dim;k++) {
         if(mask.getLabelValue(i,j,k) > 0) {
           if(getVoxelValue(i,j,k) > maximum) maximum = getVoxelValue(i,j,k);
         }
       } 
    }
  }
  return(maximum);
}

float FloatImage::imageMin(const LabelImage& mask) const
{   
  int i,j,k;
  float minimum;
  
  minimum = MAXFLOAT;
  for(i = 0;i < header.x_dim;i++) {
    for(j = 0;j < header.y_dim;j++) {
       for(k = 0;k < header.z_dim;k++) {
         if(mask.getLabelValue(i,j,k) > 0) {
           if(getVoxelValue(i,j,k) < minimum) minimum = getVoxelValue(i,j,k);
         }
       } 
    }
  }
  return(minimum);
}


//
//WARNING: this seem to be a recursive filter , is it what's expected!?
void FloatImage::knnFilter(const LabelImage& mask, int windowLenX, int windowLenY, int windowLenZ, int knn)
{
  int x,y,z,i,j,k;
  int n;
  float f,sum,voxelVal;
  char c;

  std::vector<float> window(knn);
  std::vector<float> distances(knn);

  for(x = 0;x < header.x_dim;x++) {
    for(y = 0;y < header.y_dim;y++) {
      for(z = 0;z < header.z_dim;z++) {
        if(mask.getLabelValue(x,y,z) >= 1) {
          n = 0;
          voxelVal = getVoxelValue(x,y,z);
          for(i = (-windowLenX);i < (windowLenX + 1);i++) {
            for(j = (-windowLenY);j < (windowLenX + 1);j++) {
              for(k = (-windowLenZ);k < (windowLenX + 1);k++) { 
                if(n < knn) {
                  window[n] = getSafeVoxelValue(x + i, y + j, z + k);
                  distances[n] = fabs(window[n] - voxelVal);
                  n++;
                }
                else {
                  c = maxArg(&distances[0],knn);
                  f = fabs(getSafeVoxelValue(x + i, y + j, z + k) - voxelVal);
                  if(f < distances[c]) {
                    distances[c] = f;
                    window[c] = getSafeVoxelValue(x + i, y + j, z + k);
                  }
                }
              }
            }
          }
          sum = 0.0;
          for(i = 0;i < knn;i++) {
            sum = sum + window[i];
          }
          putVoxelValue(x,y,z,sum/((float) knn));
        }
      }  
    }
  }
}

unsigned char LabelImage::labelMax(void)
{   
  int i,j,k;
  unsigned char maximum;
  
  maximum = 0;
  for(i = 0;i < header.x_dim;i++) {
    for(j = 0;j < header.y_dim;j++) {
       for(k = 0;k < header.z_dim;k++) {
         if(getLabelValue(i,j,k) > maximum) maximum = getLabelValue(i,j,k);
         
       } 
    }
  }
  return(maximum);
}

  // 3D erosion with 3 x 3 x 3 structuring element
void LabelImage::erode3D(void) 
{
  int i,j,k,x,y,z;
  bool boundaryVoxel;
  LabelImage tmp;

  tmp.copyImage(*this);
  for(i = 0;i < header.x_dim;i++) {
    for(j = 0;j < header.y_dim;j++) {
      for(k = 0;k < header.z_dim;k++) {
        if(tmp.getLabelValue(i,j,k) >= 1) {
          boundaryVoxel = false;
          for(x = -1; x < 2; x++) {
            for(y = -1; y < 2; y++) {
              for(z = -1; z < 2; z++) {
                if(tmp.getSafeLabelValue(i + x, j + y, k + z) < 1) {
                  boundaryVoxel = true;
                }
              }
            }
          }
          if(boundaryVoxel) {
            putLabelValue(i,j,k,0);
          }
        }
      }
    }
  }
}

