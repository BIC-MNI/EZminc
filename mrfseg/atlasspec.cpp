// ******************************************************
// Functions for atlas specification files
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


#include <iostream>
#include <fstream>

#include "atlasspec.h"
#include "minc_generic_image.h"

using namespace std;
//////////////////////////////////
// the atlas file should look like this:
// lines starting with # are for comments
//////////////////////////////////
// ATLAS NAME
// n numberOfLabels
// Region1name Region1filename LabelNOTPermitted1 2 3 ...
// Region2name Region2filename LabelNOTPermitted1 2 3 ...
//  ....
// Regionnname Regionnfilename
// Label1name Label1type1 Label1type2 Label1type3
// ...
// Labelnname Labelntype1 ....  
// Mrf matrix (numberofLabels + 1) x (numberOfLabels + 1)
// should be symmetric

void constructLabelList(LabelList* labels,int maxLabel,int* remThese, int howManyToRemove) 
{
  int i,j,k;
  bool addLabel;

  k = 0;
  for(i = 1;i < (maxLabel + 1);i++) {
    addLabel = true;
    for(j = 0;j < howManyToRemove;j++) {
      if(remThese[j] == i) addLabel = false;
    }
    if(addLabel) {
      labels->list[k] = i;
      k++;
    }
  }
  labels->len = k;
} 

// inline void invertLabelList(LabelList* newlabels,LabelList* labels,int maxLabel, bool inclBG)
// {}

void printLabelInfo(AtlasSpec* atlas) 
{
  int i,j;

  cout << "Number of regions:" << atlas->n << " Number of Labels:" << atlas->numberOfLabels << endl;
  cout << "Label:   Name   Pure   compos " << endl;
  for(i = 0; i < atlas->numberOfLabels;i++) {
    cout << i << "        " << atlas->labelnames[i] << "   " << atlas->labelTypes[i].pureLabel << "   ";
    if(atlas->labelTypes[i].pureLabel) cout << endl;
    else cout <<  atlas->labelTypes[i].mixed[0] << " " <<  atlas->labelTypes[i].mixed[1] << endl;   
  }
  cout << "MRF" << endl;
  
  for(i = 0; i < atlas->numberOfLabels;i++) { 
    for(j = 0; j < atlas->numberOfLabels;j++) {
      cout.width(4); 
      cout << atlas->mrfConstants[i][j];
    }
    cout << endl;
  } 
}


//////////////////////////////////
// the atlas file should look like this:
// lines starting with # are for comments
//////////////////////////////////
// ATLAS NAME
// n numberOfLabels
// Region1name Region1filename LabelNOTPermitted1 2 3 ...
// Region2name Region2filename LabelNOTPermitted1 2 3 ...
//  ....
// Regionnname Regionnfilename
// Label1name Label1type1 Label1type2 Label1type3
// ...
// Labelnname Labelntype1 ....  
// Mrf matrix (numberofLabels + 1) x (numberOfLabels + 1)
// should be symmetric

AtlasSpec::AtlasSpec()
{
}

bool AtlasSpec::readAtlasSpec(char* filename) 
{
  char str[256];
  char* pstr = str;
  int i,j;
  char c;
  int tmp[256];
  bool probLimits;

  
  ifstream ifile;
  ifile.open(filename);
  if(!ifile.is_open()) return(false);
  ifile >> pstr;
  if(!(strcmp(pstr,"p"))) {
    probLimits = true;
  }
  else { 
    probLimits = false;   
  }
  ifile >> n;
  ifile >> numberOfLabels;
 
  if(!probLimits) {
    regionLowProb.clear();
    regionUpProb.clear();
  }
  else {
    if(n > 0) {
      regionLowProb.resize(n,std::vector<float>(numberOfLabels,0.0));
      regionUpProb.resize(n,std::vector<float>(numberOfLabels,0.0));
    }
    else { 
      regionLowProb.resize(1,std::vector<float>(numberOfLabels,0.0));
      regionUpProb.resize(1,std::vector<float>(numberOfLabels,0.0));
    }
  }   

  // reserving the memory 
  if(n > 0) {
    filenames.resize(n) ;
    regionnames.resize(n) ;
    permittedLabels.resize(n);  
  }

  labelnames.resize(numberOfLabels);
  labelTypes.resize(numberOfLabels);
  mrfConstants.resize(numberOfLabels,std::vector<float>(numberOfLabels,0.0));
  

  /*for(i = 0;i < n;i++) {
    regionnames[i] = new char[MAX_NAME_LEN];
    filenames[i] = new char[MAX_NAME_LEN];
  } 

  for(i = 0;i < numberOfLabels;i++) {
    labelnames[i] = new char[MAX_NAME_LEN];
    mrfConstants[i] = new float[atlas->numberOfLabels];
   }*/

  // set the background label info
  labelnames[0]="background";
  labelTypes[0].pureLabel = true; // pure label;  

  if(probLimits && n == 0 ) {
    regionLowProb[0][0] = 0.0;
    regionUpProb[0][0] = 0.0;
    for(j = 0; j < (numberOfLabels - 1);j++) {
      ifile >> regionLowProb[0][j + 1];
      if(regionLowProb[0][j + 1] < 0.0) regionLowProb[0][j + 1] = 0.0;
      ifile >> regionUpProb[0][j + 1];
      if( regionUpProb[0][j + 1] > 1.0) regionUpProb[0][j + 1] = 1.0;
    }
  }
  // read the region info
  for(i = 0;i < n;i++) {
    ifile >> pstr;
    regionnames[i]=pstr;
    c = ifile.get();
    ifile >> pstr;
    filenames[i]=pstr; 
    if(probLimits) {
      regionLowProb[i][0] = 0.0;
      regionUpProb[i][0] = 0.0;
      for(j = 0; j < (numberOfLabels - 1);j++) {
        ifile >> regionLowProb[i][j + 1];
        if(regionLowProb[i][j + 1] < 0.0) regionLowProb[i][j + 1] = 0.0;
        ifile >> regionUpProb[i][j + 1];
        if( regionUpProb[i][j + 1] > 1.0) regionUpProb[i][j + 1] = 1.0;
      }
      constructLabelList(&(permittedLabels[i]),numberOfLabels - 1,NULL,0); 
    }
    else {
      c = ifile.get();
      j = 0;
      if(c == ' ') {     
        while(c == ' ') {
          ifile >> tmp[j];
          c = ifile.get();   
          j++;
        }
      }      
      constructLabelList(&(permittedLabels[i]),numberOfLabels - 1,tmp,j);
    }
  } 
  // read the label info
  for(i = 1;i< numberOfLabels;i++) {
    ifile >> pstr;
    labelnames[i]=pstr;
    ifile >> labelTypes[i].pureLabel;
    if(ifile.get() != '\n') {
      ifile >> labelTypes[i].mixed[0];
      ifile >> labelTypes[i].mixed[1];
    }
  }
  
  // read the mrf matrix
  for(i = 0;i < numberOfLabels;i++) {
    for(j = 0;j < numberOfLabels;j++) {
      ifile >> mrfConstants[i][j];
    }
  }
  
  ifile.close();
  return true;
}

void AtlasSpec::freeAtlas() 
{
}

AtlasSpec::~AtlasSpec()
{
  freeAtlas();
}

// reads the images in the atlas and ensures that each
// voxelsums  are equal to the unity
// returns 0 if ok
int AtlasSpec::readAtlasImages(std::vector<FloatImage>& atlasImages)
{
  int i,j,imgsize;
  int intstatus;
  int dimx,dimy,dimz;
  float probsum;
  FloatImage psum;
  bool filter = false;
  atlasImages.resize(n);
  // read images
  for(i = 0;i < n;i++) {
    std::cout<<"reading atlas image:"<<filenames[i].c_str()<<std::endl;
    intstatus = readImage(filenames[i].c_str(),atlasImages[i]);
    if(intstatus != 0) return(intstatus + i*100);
  }
  // check that the dimensions of the atlas match
  dimx = atlasImages[0].header.x_dim;
  dimy = atlasImages[0].header.y_dim; 
  dimz = atlasImages[0].header.z_dim;
  for(i = 1;i < n;i++) {
    if(dimx != atlasImages[i].header.x_dim) return(5 + i*10);
    if(dimy != atlasImages[i].header.y_dim) return(5 + i*10);
    if(dimz != atlasImages[i].header.z_dim) return(5 + i*10);
  }
  // make sure that every voxel sums to the unity in prob atlas
  imgsize = dimx*dimy*dimz;
  psum.copyImage(atlasImages[0]);

  for(i = 0;i < imgsize;i++) {
    probsum = 0.0;
    for(j = 0;j < n;j++) {
      probsum = probsum + atlasImages[j].data[i];
    } 
    psum.data[i] = probsum;
    if(fabs(probsum) > 0.0001) {
      for(j = 0;j < n;j++) {
        atlasImages[j].data[i] =  atlasImages[j].data[i]/probsum;
      }
    }       
  }
  if(filter) {
   for(i = 0;i < n;i++) {
      averageFilterZeros(atlasImages[i],psum,0.0001,1);
    }
    for(i = 0;i < imgsize;i++) {
      probsum = 0.0;
      for(j = 0;j < n;j++) {
        probsum = probsum + atlasImages[j].data[i];
      } 
   
      if(fabs(probsum) > 0.0001) {
        for(j = 0;j < n;j++) {
          atlasImages[j].data[i] =  atlasImages[j].data[i]/probsum;
        }
      }       
    } 
  }   
  return(0);
}

// takes a mask and masks the atlas (makes atlas values out of mask 0)
// in addition tries to ensure that each voxel within the mask is represented in the atlas
// mask should have value 0 for background and 1 for foreground (i.e. brain) 
// atlas images should be normalized to sum to the unity.

 
bool AtlasSpec::maskAtlas(std::vector<FloatImage>& atlasImages,const FloatImage& mask)
{
    FloatImage probmask;
    int i,j,k,l,ci,cj,ck,imgsize;
    
    // cout << "Entering the maskatlas" << endl;
    // cout << atlasImages[0]->header.x_dim << endl; 
    // cout <<  mask->header.x_dim << endl;
    // cout << atlasImages[0]->header.y_dim << " " <<  mask->header.y_dim << endl;
    // cout << atlasImages[0]->header.z_dim << " " <<  mask->header.z_dim << endl;
    // check that atlas images correspond to the mask image
   if((atlasImages[0].header.x_dim != mask.header.x_dim) ||(atlasImages[0].header.y_dim != mask.header.y_dim) 
       || (atlasImages[0].header.z_dim != mask.header.z_dim)) {
     if((atlasImages[0].header.x_dim < mask.header.x_dim) ||(atlasImages[0].header.y_dim < mask.header.y_dim) 
	|| (atlasImages[0].header.z_dim <  mask.header.z_dim)) {
       cout << "ERROR: Atlas dimensions do not match to the maskfile dimensions" << endl;
       return(false);
     }
     else {
       cout << "Warning:  Atlas dimensions do not match to the maskfile dimensions -> continuing but check the results" << endl;
     }
   }
   probmask.newImage(atlasImages[0]);
   cout << "probmask created ok" << endl;
   imgsize = (atlasImages[0].header.x_dim)*(atlasImages[0].header.y_dim)*(atlasImages[0].header.z_dim); 
   for(i = 0;i < imgsize;i++) {
     probmask.data[i] = 0;
     for(j = 0;j < n;j++) {
       probmask.data[i] = probmask.data[i] + atlasImages[j].data[i];
     }
   }
   // writeImage("tmpprob",&probmask,true);
   //cout << getVoxelValue(mask,190,32,105) << endl;
   // cout << getVoxelValue(&probmask,190,32,105) << endl;

   cout << "probmask generated" << endl;
   for(i = 0;i < mask.header.x_dim;i++ ) {
     for(j = 0;j < mask.header.y_dim;j++ ) {
       for(k = 0;k < mask.header.z_dim;k++ ) {
         // first check if the atlas spatially cover the mask
         if((mask.getVoxelValue(i,j,k) > 0.5) && (probmask.getVoxelValue(i,j,k) < 0.5) ) { 
	   //  cout << "*";
           if(probmask.findClosestNonZero(i,j,k,&ci,&cj,&ck,0.5) == false) {
             cout << "here" << i  << " " << j << " " << k<< endl;
             return(false);
           }
           float v=atlasImages[l].getVoxelValue(ci,cj,ck);
           for(l = 0;l < n;l++) {
              atlasImages[l].putVoxelValue(i,j,k,v);
           }
         }
         
         // then remove those parts of the atlas that are not necessary
         if(mask.getVoxelValue(i,j,k) < 0.5) {
           for(l = 0;l < n;l++) {
              atlasImages[l].putVoxelValue(i,j,k,0.0);
           }
         }
       }
     }
   }
   return(true);
}
