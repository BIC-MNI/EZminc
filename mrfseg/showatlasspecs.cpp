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



#include "atlasspec.h"


int main(int argc, char** argv) 
{
  AtlasSpec atlas;
  bool boolstatus;

  if(argc < 2) {
    cout << "Usage: showatlaspecs atlasfile" << endl;
    return(2);
  } 

  boolstatus = readAtlasSpec(&atlas,argv[1]);
  if(boolstatus == false) {
    cout << "Could not read atlas file " << argv[1] << endl;
    return(1);
  }  
  atlas.n = 1;
  printLabelInfo(&atlas);
  return(0);
} 

  
