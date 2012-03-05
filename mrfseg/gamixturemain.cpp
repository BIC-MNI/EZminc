// ******************************************************
// Mixture model optimization using genetic algorithms
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
// Version 1.0: The first version of the software.
// Version 1.1: Added the possibility of restarts, and speeded up the Parzen window computation. JT 17th Jul 2007.
 
#include <vector>
#include "gamixture.h"
#include "atlasspec.h"
#include "minc_generic_image.h"

int main(int argc,char** argv)
{
  FloatImage img;
  FloatImage mask;
  LabelImage mask2;
  std::vector<FloatImage> atlasImages;
  MixtureSpec mixture;
  AtlasSpec atlas; 
  PdfEstimate hatf;
  Parameters params;  
  Population pop,selPop,popRuns;

  float maxMu,minMu,minVar,maxVar,mean;
  std::vector<float> lowLimit;
  std::vector<float> upLimit;  
  //  float* fitness;

  bool terminate;
  bool boolstatus;
  bool changed_n = false;

  int intstatus,i,j,itercount,rr,n; 
  int pveLabels,pureLabels;  

  if(argc < 5) {
    cout << "Genetic algorithm for mixture model optimization" << endl;
    cout << "Usage: gamixture [-unsigned] imagefile maskfile atlasfile fmmparamfile [paramaters]" << endl;
    cout << "Optional parameters are :" << endl;
    cout << "-alpha:             parameter for blended crossover, defaults to " << DEFAULT_ALPHA << endl;
    cout << "-size:              population size, defaults to " << DEFAULT_POPSIZE << endl;
    cout << "-terminationthr:    threshold for terminating the algorithm, defaults to "<< DEFAULT_TERMINATIONTHR << endl; 
    cout << "-xoverrate:         crossover rate, defaults to " << DEFAULT_XOVERRATE << endl;
    cout << "-maxgenerations:    maximum number of generations, defaults to " << DEFAULT_MAXGENERATIONS << endl;
    cout << "-sortpop            whether to use the permutation operator, defaults to " << DEFAULT_SORTPOP << endl;
    cout << "-parzenn            number of points for Parzen estimate (for approximate ML), defaults to " << DEFAULT_PARZENN << endl;
    cout << "-parzensigma        window width parameter for the Parzen estimate, defaults to " << DEFAULT_PARZENSIGMA << endl;         
    cout << "-equalvar           whether component densities should have equal variances, defaults to " << DEFAULT_EQUALVAR << endl;  
    cout << "-restarts           the number of individual runs of the GA, defaults to " << DEFAULT_RESTARTS << endl;  
   return(1);
  }
  intstatus = readImage(argv[1],img);
  if(intstatus != 0) {
    cout << "Could not read image file " << argv[1] << " " << intstatus << endl;
    return(2);
  }  
  cout << "Image dimensions: " << img.header.x_dim << " " << img.header.y_dim << " " << img.header.z_dim << endl;
  if(!(strcmp(argv[2],"default"))) {
    mask.copyImage(img);
    mask.thresholdImage(VERY_SMALL);
    mask2.thresholdImage(img,VERY_SMALL); //!?!?!?!
  }
  else {
    intstatus = readImage(argv[2],mask);
    if(intstatus != 0) {
      cout << "Could not read image (brainmask) file " << argv[2] << " " << intstatus << endl;
      return(3);
    }  
    mask.thresholdImage(0.0001);
    mask2.thresholdImage(mask,0.5);//!?!?!
  } 
   
  boolstatus = atlas.readAtlasSpec(argv[3]);
  if(boolstatus == false) {
    cout << "Could not read atlas file " << argv[3] << endl;
    return(4);
  }   
 
  intstatus = params.parseParams(argc - 4,&(argv[4]));
  if(intstatus != 0) {  
    cout << "Incorrect parameter value input" << endl;
    return(6);
  }
  if(atlas.n > 0) {
    atlasImages.resize(atlas.n);
    intstatus = atlas.readAtlasImages(atlasImages);
    if(intstatus != 0) {
      cout << "Could not read probabilistic atlas. Error: " << intstatus << endl;
    return(6);
    }  
    cout << "The atlas files have been read" << endl;
    if(atlas.maskAtlas(atlasImages,mask) == false) {
      cout << "Could not mask atlas" << endl;
      return(7);
    }  
    cout << "The atlas has been masked" << endl;
  }
  else {
   // taking care for the case where atlas->n == 0
   // then we assume that the atlas is defined by the mask
    atlas.n = 1;
    changed_n = true;
    atlasImages.resize(atlas.n);
    boolstatus = atlasImages[0].copyImage(mask);
  } 
  mixture.allocateMixtureSpec(atlas);
  
  // compute the number of pve labels and pure labels.
  // remember that pve labels have to have indeces greater than pure labels
  pureLabels = 0;
  pveLabels = 0;
  for(i = 0;i < atlas.numberOfLabels;i++) {
    if(atlas.labelTypes[i].pureLabel == true) pureLabels++;
    else pveLabels++;     
  }
  cout << "purelabels:  " << pureLabels << " Pvelabels: " << pveLabels << endl; 
  
  
  // initializing the Parzen windows
  hatf.n = params.parzenN;
  computeX(&hatf,img,mask2,0.999);
  setSigma(&hatf,params.parzenSigma); 
  
  // Setting the limits for GA
  maxMu = hatf.x[hatf.n - 1];
  minMu = hatf.x[0];
  minVar = hatf.sigma * hatf.sigma;
  maxVar = 0.0;  // this is re-set later on  

  lowLimit.resize(3*pureLabels + pveLabels); 
  upLimit.resize(3*pureLabels + pveLabels); 
  // initialize by defults, then adapt the necessary values
  for(i = 0; i < pureLabels;i++) {
    upLimit[3*i] = maxMu;
    lowLimit[3*i] = minMu;
    upLimit[3*i + 1] = maxVar;  // this is re-set later on
    lowLimit[3*i + 1] = minVar;
    if(changed_n)  {
      upLimit[3*i + 2] = 1.0;  
    }
    else { 
      upLimit[3*i + 2] = 0.0; 
    } // this is changed afterwards
    lowLimit[3*i + 2] = 0.0;
  }
  upLimit[0] = 0.0;  // adjustment for the background label
  lowLimit[0] = 0.0;
  upLimit[2] = 0.0;
   
  // pve classes
  for(i = 0; i <pveLabels;i++) {
    if(changed_n) {
      upLimit[pureLabels*3 + i] = 1.0;
    }
    else {
      upLimit[pureLabels*3 + i] = 0.0;
    }
    lowLimit[pureLabels*3 + i] = 0.0;
  }  
  for(i = 0; i < atlas.n;i++) {
    rr = i + 1;
    cout << "Region " << rr << ":computing parzen estimate..." << endl;;  
    computeY(&hatf,img,atlasImages[i]);
    //  if(i == 5) writeEstimate("tmp.parzen",&hatf,true);
  
    cout << "Region " << rr << ":Genetic algorithm..." << endl;
    // set the region wise limits
    maxVar = 0.0;
    mean = 0.0;
    for(j = 0;j < hatf.n; j++) {
      mean = mean +  (hatf.x[1] - hatf.x[0])*(hatf.y[j])*(hatf.x[j]);
    }
    for(j = 0;j < hatf.n; j++) {
      maxVar = maxVar + (hatf.x[1] - hatf.x[0])*hatf.y[j]*(hatf.x[j] - mean)*(hatf.x[j] - mean);
    }
    maxVar = maxVar/pureLabels;
    for(j = 0; j < pureLabels;j++) {
       upLimit[3*j + 1] = maxVar;
    }
    if(!atlas.regionLowProb.empty()) {
      for(j = 1; j < pureLabels;j++) {
        lowLimit[j*3 + 2] = atlas.regionLowProb[i][j];
        upLimit[j*3 + 2] = atlas.regionUpProb[i][j];
      }
      for(j = 0; j < pveLabels;j++) {
        lowLimit[pureLabels*3 + j] = atlas.regionLowProb[i][pureLabels + j];
        upLimit[pureLabels*3 + j] = atlas.regionUpProb[i][pureLabels + j];
      }
    }
    else {
      if(!changed_n) {
        for(j = 1; j < pureLabels;j++) { // skip background
          upLimit[3*j + 2] = 0.0; 
        }
        for(j = 0; j < pveLabels;j++) {
          upLimit[pureLabels*3 + j] = 0.0;
        }
        for(j = 0;j < atlas.permittedLabels[i].len;j++) {
	  //   cout << atlas.permittedLabels[i].list[j] << " ";  
          if(atlas.permittedLabels[i].list[j] > (pureLabels -1)) {
            upLimit[3*pureLabels + atlas.permittedLabels[i].list[j] - pureLabels] = 1.0;
          }
          else {
            upLimit[atlas.permittedLabels[i].list[j]*3 + 2] = 1.0;
          }
        }
      }
    } 
    popRuns.gaInitializePopulation(params.restarts,1,pureLabels + pveLabels,pveLabels,atlas.labelTypes,lowLimit,upLimit,params.equalVar);
    for(n = 0;n < params.restarts;n++) {
      pop.gaInitializePopulation(params.size,1,pureLabels + pveLabels,pveLabels,
                           atlas.labelTypes,lowLimit,upLimit,params.equalVar);
      pop.gaSortPopulation(1);
      pop.gaEvaluate(&hatf);
      pop.gaReorder();
      selPop.copyPartialPopulation(&pop);
      terminate = false;
      itercount = 0;
      while((!terminate) && (itercount < params.maxGenerations)) {
        gaTournamentSelection(&selPop,&pop,1);
        gaBLX(&pop,&selPop,params.xoverRate,1,params.alpha,lowLimit,upLimit,params.equalVar); 
        if(params.sortPop) {
          pop.gaSortPopulation(1);
        }
        pop.gaEvaluate(&hatf);
        pop.gaReorder();
        terminate = pop.gaTerminate(params.terminationThr);
        itercount++;
      }
      if(!(params.sortPop)) pop.gaSortPopulation(1); 
      cout << "GA converged after " << itercount << " iterations." << endl;
      cout << "the KL distance is " << pop.energies[0] << "." << endl;
      // put the best invividual into popRuns
     for(j = 0;j < pureLabels;j++) {
        popRuns.gaSetMu(n,j,pop.gaGetMu(0,j));
        popRuns.gaSetSigma2(n,j,pop.gaGetSigma2(0,j));
        popRuns.gaSetProb(n,j,pop.gaGetProb(0,j));     
      }
      for(j = 0;j < pveLabels;j++) {
        popRuns.gaSetProb(n,j + pureLabels,pop.gaGetProb(0,j + pureLabels));
      }
      popRuns.energies[n] = pop.energies[0];
    }
    if( params.restarts > 1) popRuns.gaReorder(); 
    // convert the best indivual to mixtureSpec
    for(j = 0;j < pureLabels;j++) {
      mixture.putMu(i,j,popRuns.gaGetMu(0,j));
      mixture.putSigma2(i,j,popRuns.gaGetSigma2(0,j));
      mixture.putProb(i,j,popRuns.gaGetProb(0,j));     
    }
    for(j = 0;j < pveLabels;j++) {
      mixture.putProb(i,j + pureLabels,popRuns.gaGetProb(0,j + pureLabels));
    }
  }
  mixture.printMixture();
  intstatus = writeMixtureParameters(argv[4],mixture, true);
  if(intstatus != 0) {
    cout << "Error in writing the mixture parameters file" << argv[4] << endl;
    return(7);
  }
  return(0); 
}
