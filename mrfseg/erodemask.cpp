// 
// ******************************************************
// Morpholgical erosion
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


#include "analyze.h"

int main(int argc, char** argv) 
{
  AnalyzeImage img;
  AnalyzeLabelImage mask;
  int intstatus,i;
  bool signedData = true;

  if(argc > 1) {
    if(!(strcmp(argv[1],"-unsigned"))) {
      signedData = false;
      argc--;
      for(i = 1;i < argc;i++) {
        argv[i] = argv[i + 1];
      }
    }
  }
  if(argc < 3) {
    cout << "Usage: erodemask [-unsigned] imagefilename eroded_mask_filename" << endl; 
    return(1);
  }

  cout << argv[1] << endl;
  intstatus = readImage(argv[1],&img,signedData);
  if(intstatus != 0) {
    cout << "Could not read image file " << argv[1] << " " << intstatus << endl;
    return(2);
  }  
  binaryMask(&img,&mask,0.0001);
  erode3D(&mask);
  intstatus = writeLabelImage(argv[2],&mask,true);
  if(intstatus != 0) {
    cout << "Could not write labeled image. Error: " << intstatus << endl;
    return(3);
  } 
  
  return(0);
} 
