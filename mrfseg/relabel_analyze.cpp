// 
// ******************************************************
// Relabeling a labeled Analyze image
// (C) Jussi Tohka, 2005 
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
// *****************************************************

#include "generic_image.h"
#include "minc_generic_image.h"
#include <iostream>
using namespace std;

int main(int argc,char** argv)
{
  
  //todo: what is going on here!? 
  // shouldn't it be LabelImage s?
  LabelImage imgin,imgout;
  int x,y,z,labelsToChange,i;
  int intstatus;  
  int inlabel,outlabel;
  bool clobber = true;
  bool status;
  div_t divtmp;

  if(argc < 5) {
    cout << "Usage: relabel_analyze inputfile outputfile inlabel1 outlabel1" << endl;
    return(1);
  }
  if(readImage(argv[1],imgin)) {
    cerr << "Could not read the file " << argv[1] << endl;
    return(2);
  }
  status = imgout.copyImage(imgin);
  divtmp = div(argc - 3,2);
  if(divtmp.rem == 0) {
    labelsToChange = divtmp.quot;
  }
  else {
    cout << " wrong number of inputs" << endl;
    return(3);
  }
  for(i = 0;i < labelsToChange;i++) {
    inlabel = std::stoi(argv[2*i + 3]);
    outlabel = std::stoi(argv[2*i + 4]);
    cout << "Changing" << inlabel << "->" << outlabel << endl;
    for(int j=0;j<imgin.data.size();j++)
    {
      if(imgin.data[j]==inlabel)
            imgout.data[j]=outlabel;
    }
  }
  if(writeImage(argv[2],imgout)) {
    cerr << "Could not write the file " << argv[2] << endl;
    return(3);
  }
  return(0);
}
