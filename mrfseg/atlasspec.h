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


#ifndef ATLASSPEC_H
#define ATLASSPEC_H

#include <string>
#include <vector>

#include "generic_image.h"

#define MAX_NAME_LEN 256
#define MAX_LABELS 20 
#define PURELABEL  0
#define PVELABEL 1


struct LabelType {
  bool pureLabel;
  int mixed[2];  
};

struct LabelList {
  int len;
  int list[MAX_LABELS];
};

class AtlasSpec {
  public: //for now
  int n;                      // Gives the number of different fuzzy regions
  int numberOfLabels;         // Gives the number of possible labels, background included 
   
  std::vector<std::string> labelnames;          // Gives the names of labels; 0 is reserved for background
  std::vector<std::string> filenames;           // Gives the names of the files containg fuzzy masks
  std::vector<std::string> regionnames;         // Gives the description of each region (optional)
  std::vector<LabelList>   permittedLabels;     // Labels permitted for each region 
                       
  std::vector<std::vector<float> > mrfConstants;    // pairwise interactions in the mrf
  std::vector<LabelType> labelTypes;                // types of labels (pve or gaussian)
  std::vector<std::vector<float> > regionLowProb;
  std::vector<std::vector<float> > regionUpProb;                  
  
  public:
  AtlasSpec();
  ~AtlasSpec();
  
  bool readAtlasSpec(char* filename);
  void freeAtlas();
  int  readAtlasImages(std::vector<FloatImage>& atlasImages); 
  bool maskAtlas(std::vector<FloatImage>& atlasImages,const FloatImage& mask);
};

void constructLabelList(LabelList* labels,int maxLabel,int* remThese, int howManyToRemove) ;
void printLabelInfo(AtlasSpec* atlas) ;

#endif

