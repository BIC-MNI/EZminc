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


#include <algorithm>
#include "gamixture.h"
#include <ctime>

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


// This function parses upperLimit and lowerLimit from the command line 
// input to the program. 

int Parameters::parseParams(int n,char** arguments) 
{
  // Initialize with the default values
  int i;  

  alpha = DEFAULT_ALPHA;
  size  = DEFAULT_POPSIZE;
  terminationThr = DEFAULT_TERMINATIONTHR;
  xoverRate = DEFAULT_XOVERRATE;
  maxGenerations = DEFAULT_MAXGENERATIONS;
  sortPop   = DEFAULT_SORTPOP;
  parzenN  = DEFAULT_PARZENN;
  parzenSigma = DEFAULT_PARZENSIGMA;
  equalVar = DEFAULT_EQUALVAR;
  restarts = DEFAULT_RESTARTS;

  if( ((n - 1) % 2) > 0 ) return(1);
  for(i = 1;i < (n - 1)/2 + 1;i++) {
    if(!strcmp(arguments[2*i - 1],"-alpha"))
      alpha = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-size"))
      size = std::stoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-terminationthr"))
      terminationThr = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-xoverrate"))
      xoverRate = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-maxgen"))
      maxGenerations = std::stoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-sortpop"))
      sortPop = std::stoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-parzenn"))
      parzenN = std::stoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-parzensigma"))
      parzenSigma = std::stoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-equalvar"))
      equalVar = std::stoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-restarts"))  // Version 1.1
      restarts = std::stoi(arguments[2*i]);
  }
  return(0);
}

void displayParameterHelp()
{

}

int Population::gaInitializePopulation(int _size, int _dim, int _numberOfLabels,int _numberOfPveLabels,
                           const std::vector<LabelType>& _labelTypes, 
                           const std::vector<float>& _lowLimit, 
                           const std::vector<float>& _upLimit, bool equalVar)
{
  int i,j;
  int mixtureSize;
  int bgVar;
  float tmpRand;  
  float probsum;

  size = _size;
  dim =  _dim;
  numberOfLabels = _numberOfLabels;
  numberOfPveLabels = _numberOfPveLabels;
  labelTypes = _labelTypes;
  lowLimit = _lowLimit;
  upLimit = _upLimit;

  
  // compute the size of each mixture
  if (dim == 1) { 
    mixtureSize = 2*dim*(numberOfLabels - numberOfPveLabels) + numberOfLabels;
  }
  else {
    mixtureSize = (2*dim + dim*dim)*(numberOfLabels - numberOfPveLabels) + numberOfLabels;
  }
  // compute the number of background variables
  if (dim == 1) {
    bgVar = 3;
  }
  else {
    bgVar = 2*dim + 1 + dim*dim;
  } 
  // Allocate
  mixtures.resize(size,std::vector<float>(mixtureSize,0.0));
  energies.resize(size,0.0);
  
  srand(time(0));
  for(i = 0;i < size;i++) {
    gaSetMu(i,0,0.0);    // background
    gaSetProb(i,0,0.0);  // background
    for(j = bgVar; j < mixtureSize;j++) {
      tmpRand = (float) rand() / RAND_MAX;
      mixtures[i][j] = tmpRand*(_upLimit[j] - _lowLimit[j]) + _lowLimit[j];
    }
    gaSetSigma2(i,0,gaGetSigma2(i,1));  // background variance
    if(equalVar) {
      for(j = 2;j < (numberOfLabels - numberOfPveLabels);j++) {
        gaSetSigma2(i,j,gaGetSigma2(i,1));
      }
    }
    probsum = 0.0;
    for(j = 1;j < numberOfLabels;j++) {
      probsum = probsum + gaGetProb(i,j);
    }
    for(j = 1;j < numberOfLabels;j++) {
      gaSetProb(i,j,gaGetProb(i,j) / probsum);
    }
  }
  return(0);
}


// computes the Kullback Leibler divergence between the mixture and the 
// given data density
// ASSUMES THAT X-AXIS OF THE DENSITY IS EQUALLY SPACED

void Population::gaEvaluate(PdfEstimate* hatf)
{
  int i,j,k;
  float interval;
  float klAddEnergy;
  float mixtureVal;
  float tmp;
  int pveIntervals = 10;

  klAddEnergy = 0.0;
  interval = hatf->x[1] - hatf->x[0];
  for(i = 0;i < hatf->n ; i++) {
    if ( hatf->y[i] > 0 ) {
      klAddEnergy = klAddEnergy + hatf->y[i]*log( hatf->y[i] );
    }
  }
  klAddEnergy = interval*klAddEnergy;
  
  for(i = 0;i < size;i++) {
    energies[i] = 0.0;
    for(j = 0;j < hatf->n;j++) {
      mixtureVal = 0.0;
      for(k = 1;k < (numberOfLabels - numberOfPveLabels);k++) {
        mixtureVal = mixtureVal + gaGetProb(i,k)*computeGaussianLikelihood(hatf->x[j],gaGetMu(i,k),gaGetSigma2(i,k));
       
      }
      
      for(k = (numberOfLabels - numberOfPveLabels);k < numberOfLabels;k++) {
        mixtureVal = mixtureVal + gaGetProb(i,k)*computeMarginalizedLikelihood(hatf->x[j],gaGetMu(i,labelTypes[k].mixed[0]), 
                                                                                   gaGetMu(i,labelTypes[k].mixed[1]), 
                                                                                   gaGetSigma2(i,labelTypes[k].mixed[0]), 
                                                                                   gaGetSigma2(i,labelTypes[k].mixed[1]),10 ,DEFAULT_PVESTART, DEFAULT_PVEEND);
      }
      if(mixtureVal == 0) {
        mixtureVal = 1/(interval*1000000);
      }
      energies[i] = energies[i] + hatf->y[j]* log(mixtureVal);
    }
    energies[i] = energies[i]*interval; 
    energies[i] = klAddEnergy - energies[i];
  }
  
}


// selects pop->size - elitism from the population by means of the tournament selection
// The tournament size is 2 in the current implementation.
// Assumes that the population has been sorted.
 void gaTournamentSelection(Population* selPop,Population* pop,int elitism)
{
  int i,j,mixtureSize,r,r1,r2;
 

  // compute the size of each mixture
  if (pop->dim == 1) { 
    mixtureSize = 2*(pop->dim)*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  else {
    mixtureSize = (2*(pop->dim) + (pop->dim)*(pop->dim))*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  // copyPartialPopulation(pop,selPop);
  srand(time(0));
  for(i = 0;i < elitism;i++) {
    for(j = 0;j < mixtureSize;j++) {
      selPop->mixtures[i][j] = pop->mixtures[i][j];  
    }                                                
    selPop->energies[i] = pop->energies[i];
  }    
  for(i = elitism;i < pop->size;i++) {
    r1 = rand() % pop->size;
    r2 = rand() % pop->size;
    r = MIN(r1,r2);
    for(j = 0;j < mixtureSize;j++) {
      selPop->mixtures[i][j] = pop->mixtures[r][j];  // this works because population is 
    }                                                // ordered according to fitness
    selPop->energies[i] = pop->energies[r];
  } 
}

void Population::gaReorder(void)
{
  int i,j,mixtureSize;
  SortStruct* ss;
  
  std::vector<std::vector<float> > mixTmp;
 
  if (dim == 1) { 
    mixtureSize = 2*(dim)*(numberOfLabels - numberOfPveLabels) + numberOfLabels;
  }
  else {
    mixtureSize = (2*(dim) + (dim)*(dim))*(numberOfLabels - numberOfPveLabels) + numberOfLabels;
  }
 
  ss = new SortStruct [size];
  for(i = 0;i < size;i++) {
    ss[i].energy = energies[i];
    ss[i].index = i;
  }
  mixTmp.resize(size,std::vector<float>(mixtureSize,0.0)); //shouldn't even do this...
  
  //todo:replace this with std::sort
  qsort(ss, size, sizeof(SortStruct), cmpStruct);
  mixTmp=mixtures;
  
  for(i = 0;i < size;i++) {
    energies[i] = ss[i].energy;
    for(j = 0;j < mixtureSize;j++) {
      mixtures[i][j] = mixTmp[ss[i].index][j];
    }
  }
  delete[] ss;
}

// Implements the blended crossover
// Assumes that population pop is sorted and 
// that the selPop is generated by the tournament selection
void gaBLX(Population* pop,Population* selPop, 
           float xoverRate, int elitism, float alpha, 
           const std::vector<float>& lowLimit, const std::vector<float>& upLimit, bool equalVar)
{
  int i,j;
  int leaveAlone;
  int mixtureSize;
  int bgVar;
  int par1,par2;
  float b,probsum;

  // first skip the indivuals based on the elistism
  // i.e elistism first (that is fittest) indivuals survive 
  // automatically and unchanged to the next generation

  // then if xoverRate < 1.0 copy necessary number of indivuals from selPop
  // to pop. this is not the typical way to implement crossover but it nevertheless 
  // works because the selPop is in random order
  if (pop->dim == 1) { 
    mixtureSize = 2*(pop->dim)*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  else {
    mixtureSize = (2*(pop->dim) + (pop->dim)*(pop->dim))*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  if (pop->dim == 1) {
    bgVar = 3;
  }
  else {
    bgVar = 2*pop->dim + 1 + (pop->dim)*(pop->dim);
  } 
  leaveAlone = (int) floor( (1 - xoverRate) * (pop->size));
  for(i = elitism;i < elitism + leaveAlone;i++) {
    for(j = 0;j < mixtureSize;j++) {
      pop->mixtures[i][j] = selPop->mixtures[i][j];
    }
  }
  srand(time(0));
  for(i = elitism + leaveAlone;i < pop->size;i++) {
    par1 = rand() % pop->size;
    par2 = rand() % pop->size;
    for(j = bgVar;j < mixtureSize;j++) {
      b = (float) rand()/RAND_MAX;
      b = (1 + 2*alpha)*b - alpha;
      if((b < (-0.5)) || (b > 1.5)) cout << "b" << b << endl; 
      pop->mixtures[i][j] = b*(selPop->mixtures[par1][j]) + (1 - b)*(selPop->mixtures[par2][j]);
      pop->mixtures[i][j] = MIN(pop->mixtures[i][j],upLimit[j]);
      pop->mixtures[i][j] = MAX(pop->mixtures[i][j],lowLimit[j]);
    }
    probsum = 0.0;
    for(j = 1;j < pop->numberOfLabels;j++) {
      probsum = probsum + pop->gaGetProb(i,j);
    }
    for(j = 1;j < pop->numberOfLabels;j++) {
      pop->gaSetProb(i,j,pop->gaGetProb(i,j) / probsum);
    }
    // take care of the background variance
    pop->mixtures[i][1] = pop->mixtures[i][4];
    if(equalVar) {
      for(j = 2;j < (pop->numberOfLabels - pop->numberOfPveLabels);j++) {
        pop->gaSetSigma2(i,j,pop->gaGetSigma2(i,1));
      }
    }
  }
}

// Implements the sort operator. Uses the trivial Bubble Sort to achieve its goal.
// The Bubble Sort should be handy because the number of floats that need to be sorted is small.
void Population::gaSortPopulation(int sortDim)
{
  int i;
  int componentSize,pureLabels;
  
  pureLabels = numberOfLabels - numberOfPveLabels;
   if (dim == 1) { 
    componentSize = 3;
  }
  else {
    componentSize = 2*(dim) + (dim)*(dim) + 1;
  }  
  for(i = 0;i < size;i++) {
    //todo: change this to std::sort
    bubbleSort(&mixtures[i][0],pureLabels,componentSize);
  }
}

// whether to terminate
// Assumes sorted population
// return 'true' if termination condition is fulfilled

bool Population::gaTerminate(float thr) 
{
  float mean,best;
  int i;
  
  mean = 0.0;
  for(i = 0;i < size;i++) {
    mean = mean + energies[i];
  }
  mean = mean / size;
  best = energies[0];
  return( (mean - best) < thr);
}




