// ******************************************************
// MRFSEG main function
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


#include "generic_image.h"
#include "minc_generic_image.h"
#include "atlasspec.h"
#include "nhmrf.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv) 
{
  FloatImage img;
  FloatImage mask;
  LabelImage labelImg;
  LabelImage pveLabelImg;
  std::vector<FloatImage> atlasImages;
  std::vector<FloatImage> labelLikelihoods;
  MixtureSpec mixture;
  AtlasSpec atlas;
  
  int intstatus,i,itercount; 
  
  bool boolstatus,clobber,markov;

  //  bool em;

  float beta1,beta2;

  // the call should be :        mrfseg imagefilename brainmask_filename
  //                             atlas_def_filename mixture_params_filename 
  //                             labelfilename [pvelabelfilename] [beta1] [beta2] [markov]
  if(argc < 6) {
    cout <<  "Usage: mrfseg img brainmask atlas_def mixture_params labelimg [pvelabelimg] [beta1] [beta2]" << endl;
    cout <<  "If you say 'default' for brainmask, it is assumed that all non zero intensities in" << endl;
    cout <<  "input image are within the brain mask" << endl;  
    return(1); 
  }

  intstatus = readImage(argv[1],img);
  if(intstatus != 0) {
    cout << "Could not read image file " << argv[1] << " " << intstatus << endl;
    return(2);
  }  
  if(!(strcmp(argv[2],"default"))) {
    mask.copyImage(img);
    mask.thresholdImage(VERY_SMALL);
  }
  else {
    if(intstatus=readImage(argv[2],mask)) {
      cerr << "Could not read image (brainmask) file " << argv[2] << " " << intstatus << endl;
      return(3);
    }  
    mask.thresholdImage(0.0001);
  }
  if(!atlas.readAtlasSpec(argv[3])) {
    cerr << "Could not read atlas file " << argv[3] << endl;
    return(4);
  }  

  if(intstatus= readMixtureParameters(argv[4],atlas,mixture)) {
    cout << "Could not read mixture file " << argv[4] << " " << intstatus << endl;
    return(5);
  }  

  if(argc > 7) {
    beta1 = atof(argv[7]);
  }
  else {
    beta1 = 0.0;
  }

  if(argc > 8 ) {
    beta2 = atof(argv[8]);
  }
  else {
    beta2 = 0.1;
  }
  if(argc > 9) {
    markov = !(strcmp(argv[9],"markov"));
    cout << "M " << markov << endl;
  }
  else {
   markov  = false;
  }

 
 
 // then we assume that the atlas is defined by the mask
  atlas.n = 1;
  mixture.patlas->n = 1;
   
  atlasImages.resize(atlas.n);
  boolstatus = atlasImages[0].copyImage(mask);
    

  if(!(pveLabelImg.newImage(img))) {
    cout << "Failed to create pveLabel image" << endl;
    return(8);
  }  

  if(!(labelImg.newImage(img))) {
    cout << "Failed to create Label image" << endl;
    return(8);
  }  

  cout << "Computing the ML classification" << endl; 
  labelLikelihoods.resize(atlas.numberOfLabels - 1);
  for(i = 0;i < (atlas.numberOfLabels - 1);i++) {
    if(!(labelLikelihoods[i].newImage(img))) {
      cout << "Failed to create likelihood image" << endl;
      return(8);
    }  
  }  
  // maximum likelihood labeling.
  intstatus = computeVoxelLikelihood(mixture,img,mask,atlasImages,labelLikelihoods);
  if(intstatus != 0) {
    cout << "something wrong with the ML classification " << intstatus << endl;
    return(10);
  }
  printLabelInfo(&atlas);
  cout << "Iterated Conditional modes..." << endl;
  if(markov) {
    itercount = computeMRF(pveLabelImg, mixture,mask,labelLikelihoods,atlasImages, beta1, beta2,50,true);
  }
  else {
    itercount = computeGibbs(pveLabelImg,mixture,mask,labelLikelihoods,atlasImages, beta1, beta2,50,true);
  }
  cout << "Iterations: " << itercount << endl;
  
  //todo: delete all the labelLikelihoods?
  /*
  for(i = 0;i < (atlas.numberOfLabels - 1);i++) {
    freeImage(labelLikelihoods[i]);
  }*/
  intstatus = convertPVElabels(labelImg,pveLabelImg,img,atlasImages,mixture);
  if(intstatus != 0) {
    cout << "Conversion to pure labels did not succeed" << endl;
    return(11);
  }  
  //  cout << argv[5] << endl;
  intstatus = writeImage(argv[5],labelImg);
  if(intstatus != 0) {
    cout << "Could not write labeled image " << argv[5] << ".  Error: " << intstatus << endl;
    return(12);
  }  
  if(argc > 6) {
    intstatus = writeImage(argv[6],pveLabelImg);
    if(intstatus != 0) {
      cout << "Could not write pve labeled  " << argv[6] << ".  Error: " << intstatus << endl;
      return(13);
    }
  }
  //  cout << "here!" << endl; 
  return(0);
  // think about freeing the rest of the images as well
  // right...
}


