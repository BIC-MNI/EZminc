// ******************************************************
// Functions for mixture model optimization using genetic algorithms
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


// -------------------------------------
// At the moment this module cannot handle correlation matrices
// At the moment only for 1-d distributions

#ifndef GAMIXTURE_H
#define GAMIXTURE_H

#include <vector>
#include <iostream>
#include "generic_image.h"
#include "atlasspec.h"
#include "nhmrf.h"
#include "parzen.h"

using namespace std;


// Default parameter values
#define DEFAULT_ALPHA           0.5
#define DEFAULT_POPSIZE         100
#define DEFAULT_XOVERRATE         1
#define DEFAULT_MAXGENERATIONS  500      
#define DEFAULT_TERMINATIONTHR  0.0005
#define DEFAULT_SORTPOP         1
#define DEFAULT_PARZENN         101
#define DEFAULT_PARZENSIGMA     1
#define DEFAULT_EQUALVAR        0
#define DEFAULT_RESTARTS        1


struct SortStruct {  // I know this is clumsy, but couldn't come up with better alternative in ten minutes or so.... 
  float energy;
  int index;
};

class Population {
public: //for now  
  int size;
  int dim;            //at the moment this must be equal to 1
  int numberOfLabels; // label 0 is reserved for background, this should be counted
  int numberOfPveLabels;
  std::vector<LabelType> labelTypes; // first labels should be non PVE and then the PVE labels follow
  std::vector<float>     lowLimit;       // gives the lower limit on the probability of the label
  std::vector<float>     upLimit;        // gives the upper limit on the probability of the label
  std::vector<std::vector<float> > mixtures; 
  std::vector<float> energies; 
public:
  ~Population()
  {}
  // this is only 1-d version
  float gaGetProb(int indNumber, int classNumber) 
  { 
    int num = 3;  
    if(classNumber > (numberOfLabels - numberOfPveLabels - 1)) {
      return(mixtures[indNumber][(numberOfLabels - numberOfPveLabels)*num 
              + classNumber - (numberOfLabels - numberOfPveLabels)]);
    }
    else {
      return mixtures[indNumber][classNumber*num + num - 1];
    }
  }

  // this is only 1-d version
  float gaGetMu(int indNumber, int classNumber) 
  {
    int num = 3;
    return mixtures[indNumber][classNumber*num ];
  }

  // this is only 1-d version
  float gaGetSigma2(int indNumber, int classNumber) 
  {
    int num = 3;
    return mixtures[indNumber][classNumber*num + 1];
  }
  // this is only 1-d version
  float gaSetProb(int indNumber, int classNumber, float val) 
  { 
    int num = 3;  
    if(classNumber > (numberOfLabels - numberOfPveLabels - 1)) {
      mixtures[indNumber][(numberOfLabels - numberOfPveLabels)*num + classNumber
          - (numberOfLabels - numberOfPveLabels)] = val;
    }
    else {
      mixtures[indNumber][classNumber*num + num - 1] = val;
    }
  }

  // this is only 1-d version
  float gaSetMu(int indNumber, int classNumber,float val) 
  {
    int num = 3;
    mixtures[indNumber][classNumber*num] = val;
  }

  // this is only 1-d version
  float gaSetSigma2(int indNumber, int classNumber, float val) 
  {
    int num = 3;
    mixtures[indNumber][classNumber*num + 1] = val;
  }

  // copies "essential" parts of the population
  void copyPartialPopulation(Population* sourcePop)
  {
    int i,j,mixtureSize;

    size = sourcePop->size;
    dim  = sourcePop->dim;
    numberOfLabels    = sourcePop->numberOfLabels;
    numberOfPveLabels = sourcePop->numberOfPveLabels;
    labelTypes.clear();
    upLimit.clear();
    lowLimit.clear();

    // compute the size of each mixture
    /*
    if (sourcePop->dim == 1) { 
      mixtureSize = 2*(sourcePop->dim)*(sourcePop->numberOfLabels - sourcePop->numberOfPveLabels) + sourcePop->numberOfLabels;
    }
    else {
      mixtureSize = (2*(sourcePop->dim) + (sourcePop->dim)*(sourcePop->dim))
                      *(sourcePop->numberOfLabels - sourcePop->numberOfPveLabels) + sourcePop->numberOfLabels;
    }
    // Allocate
    mixtures.resize(sourcePop->size,std::vector<float>(sourcePop->mixtureSize,0.0));
    energies.resize(sourcePop->size);
    energies=sourcePop->energies;
    for(i = 0;i < sourcePop->size;i++) {
      for(j = 0;j < mixtureSize;j++) {
        mixtures[i][j] = sourcePop->mixtures[i][j];
      }
    }*/
    energies=sourcePop->energies;
    mixtures=sourcePop->mixtures;
    
  }
  
  int gaInitializePopulation(int size, int dim, int numberOfLabels,int numberOfPveLabels,
                             const std::vector<LabelType>& labelTypes, const std::vector<float>& lowLimit, 
                             const std::vector<float>& upLimit, bool equalVar = false);
  void gaEvaluate(PdfEstimate* hatf);
  void gaReorder(void); // Re-order individuals in the population according the fitness
  void gaSortPopulation(int sortDim);
  bool gaTerminate(float thr);
};

//todo deside where these should go
void gaTournamentSelection(Population* selPop, Population* pop,int elitism);
void gaBLX(Population* pop,Population* selPop, float xoverRate, int elitism, float alpha, 
           const std::vector<float>& lowLimit, const std::vector<float>& upLimit, bool equalVar = false);

struct Parameters {
  float alpha;
  int size;
  float terminationThr;
  float xoverRate; 
  int maxGenerations;
  int sortPop;
  int parzenN;  
  float parzenSigma;
  bool equalVar;
  int restarts;
  public:
  int parseParams(int n,char** arguments); 
};

inline int cmpStruct(const void* a, const void* b)
{
  if( (((SortStruct*) a)->energy) > (((SortStruct*) b)->energy))
    return(1);
  else
    return(-1);
}

//TODO: replace it with sort from algorythms
inline void bubbleSort(float* individual , int nofMixtures, int componentSize)
{
  int i,j,k;
  float tmp;  

  for(i = 0;i < (nofMixtures - 2);i++) {
    for(j = 1; j < (nofMixtures - 1 - i);j++) {
      if(individual[(j + 1)*componentSize] < individual[j*componentSize]) {
        for(k = 0;k < componentSize;k++) {
          tmp = individual[j*componentSize + k];
          individual[j*componentSize + k] = individual[(j + 1)*componentSize + k];
          individual[(j + 1)*componentSize + k] = tmp;
          }  
      }
    }
  }
}

void displayParameterHelp();

#endif
