// ******************************************************
// Functions for SVPA based non homogeneous MRFs 
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
#include "nhmrf.h"
#include <iostream>
#include <fstream>

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

using namespace std;

// reads the parameters for the mixture density fron a file.
// the exact format of the file is decided by the given MixtureSpec,
// more precisely AtlasSpec attached to it.
// the format is 
// region 1 params 
// region 2 params
//  ....
// region n params

// the classes that are not allowed for some region should have 
// the probability parameter value 0.0 for these regions.
void MixtureSpec::printMixture() 
{
  int i,j;
  for(i = 0;i < patlas->n;i++) {
    if(i > 1) { 
      cout << patlas->regionnames[i] << endl;
    }
    for(j = 0;j < patlas->numberOfLabels;j++) {
      if(patlas->labelTypes[j].pureLabel)
        cout << getMu(i,j) << " " << getSigma2(i,j) << " " << getProb(i,j) << endl;
      else
        cout << getProb(i,j) << endl;
    }
    cout << "----------------------------------------" << endl;
  }
}

int readMixtureParameters(char* filename, AtlasSpec& atlas, MixtureSpec& mixture)
{
  int r,l;
  float val;
  bool changed_n = false;

  ifstream ifile;
  ifile.open(filename);
  if(!ifile.is_open()) return(1); 
    
  mixture.patlas = &atlas;
  if(mixture.patlas->n == 0) {
    mixture.patlas->n = 1;
    changed_n = true;
  }

  // reserving the memory
  mixture.allocateMixtureSpec(*mixture.patlas);
  
  for(r = 0;r < mixture.patlas->n;r++) {
    for(l = 1;l < mixture.patlas->numberOfLabels;l++) {
      if(mixture.patlas->labelTypes[l].pureLabel) {
        ifile >> val;
        mixture.putMu(r,l,val);       
        ifile >> val;
        mixture.putSigma2(r,l,val);
        ifile >> val;
        mixture.putProb(r,l,val);
      }
      else {
        ifile >> val;
        mixture.putProb(r,l,val);
      }
    }
  }
  // place values for the background class
  for(r = 0;r < mixture.patlas->n;r++) {
    mixture.putMu(r,0,0.0); 
    mixture.putSigma2(r,0,mixture.getSigma2(r,1));
    mixture.putProb(r,0,0.0);
  }
  
  if(changed_n) mixture.patlas->n = 0;
  ifile.close();
  return(0);
}

int writeMixtureParameters(char* filename,MixtureSpec& mixture,bool overwrite)
{
  int r,l;
  float val;
  bool changed_n = false;

  if(! overwrite) {
    ifstream ifile;
    ifile.open(filename);
    if(ifile.is_open()) {
      ifile.close();
      return(1); 
    }
  }
  ofstream ofile;
  ofile.open(filename);
  if(!ofile.is_open()) return(2);
   if(mixture.patlas->n == 0) {
    mixture.patlas->n = 1;
    changed_n = true;
  }
  for(r = 0;r < mixture.patlas->n;r++) {
    for(l = 1;l < mixture.patlas->numberOfLabels;l++) {
      if(mixture.patlas->labelTypes[l].pureLabel) {
        ofile << mixture.getMu(r,l);
        ofile << ' ';       
        ofile << mixture.getSigma2(r,l);
        ofile << ' ';
        ofile << mixture.getProb(r,l);
        ofile << ' ';
      }
      else {
        ofile << mixture.getProb(r,l);
        ofile << ' ';
      }
    }
    ofile << endl;
  }  
  if(changed_n) mixture.patlas->n = 0;
  ofile.close();
  return(0);
}

// Computes the value of the likelihood term of an intensity value in an image for each label
// AnalyzeImage(s) pointed by labelLikelihoods should be allocated.  
// Atlas images should be in the same space than the image
// The mask is assumed to be thresholded
// A value 0 is returned if everything is ok.

int computeVoxelLikelihood(MixtureSpec& mixture,const FloatImage& img, 
    const FloatImage& mask,const std::vector<FloatImage>& atlasImages,
    std::vector<FloatImage>& labelLikelihoods)
{
  int i,j,k,imgsize;
  int x,y,z;
  int dimx,dimy,dimz;
  int label1,label2;
  float lvalue,tmpval;

  // First test that the mask size and image size and atlas size match

  //  printMixture(mixture);

  dimx = img.header.x_dim;
  dimy = img.header.y_dim; 
  dimz = img.header.z_dim;

  if(( dimx != mask.header.x_dim ) || ( dimy != mask.header.y_dim ) 
                                    || ( dimz != mask.header.z_dim ))
    return(2);

  if(( dimx != atlasImages[0].header.x_dim ) || ( dimy != atlasImages[0].header.y_dim ) 
     || ( dimz != atlasImages[0].header.z_dim )) {
    if(( dimx > atlasImages[0].header.x_dim ) || ( dimy > atlasImages[0].header.y_dim ) 
       || ( dimz > atlasImages[0].header.z_dim )) {
      return(1);
    }
    else {
      cout << "Warning: image dimensions do not match to the atlas dimensions->continuing but check the results" << endl;
    }
  }
  imgsize = dimx*dimy*dimz;
  
  // remember: the label 0 is reserved for background 
  if((mixture.patlas->n) == 0) {
    if(mask.data[i] > 0.5) {
      for(k = 1;k < (mixture.patlas->numberOfLabels);k++) {
        labelLikelihoods[k - 1].data[i] = 0.0;
      }
      for(k = 1;k < (mixture.patlas->numberOfLabels);k++) {
        if(fabs(mixture.getProb(j,k)) > 0.0001) {
          if(mixture.patlas->labelTypes[k].pureLabel) {
            lvalue = computeGaussianLikelihood(img.data[i],mixture.getMu(j,k),mixture.getSigma2(j,k));
          }
          else {
            label1 = mixture.patlas->labelTypes[k].mixed[0];
            label2 = mixture.patlas->labelTypes[k].mixed[1];
            lvalue = computeMarginalizedLikelihood(img.data[i],mixture.getMu(j,label1),
                                                       mixture.getMu(j,label2),
                                                       mixture.getSigma2(j,label1),
						       mixture.getSigma2(j,label2),10,0.0,1.0);
          }
        }
        else {
          lvalue = 0.0;
        }
        labelLikelihoods[k - 1].data[i] = labelLikelihoods[k - 1].data[i] + lvalue;      
      }
    }
  }
  else {
    //  for(i = 0;i < imgsize;i++) {
    for(x = 0;x < dimx;x++) {
      for(y = 0;y < dimy;y++) {
        for(z = 0;z < dimz;z++) {
          if(mask.getVoxelValue(x,y,z) > 0.5) {
            for(k = 1;k < (mixture.patlas->numberOfLabels);k++) {
              labelLikelihoods[k - 1].putVoxelValue(x,y,z,0.0);
            }
            for(j = 0; j < mixture.patlas->n; j++) {
              if(atlasImages[j].getVoxelValue(x,y,z) > 0.000001) {
                for(k = 1;k < (mixture.patlas->numberOfLabels);k++) {
                  if(fabs(mixture.getProb(j,k)) > 0.0001) {
                    if(mixture.patlas->labelTypes[k].pureLabel) {
                      lvalue = computeGaussianLikelihood(img.getVoxelValue(x,y,z),mixture.getMu(j,k),
                                                   mixture.getSigma2(j,k));
                    }
                    else {
                      label1 = mixture.patlas->labelTypes[k].mixed[0];
                      label2 = mixture.patlas->labelTypes[k].mixed[1];
                      lvalue = computeMarginalizedLikelihood(img.getVoxelValue(x,y,z),
                                                       mixture.getMu(j,label1),
                                                       mixture.getMu(j,label2),
                                                       mixture.getSigma2(j,label1),
                      mixture.getSigma2(j,label2),10,DEFAULT_PVESTART,DEFAULT_PVEEND);
                    }
                  } else {
                    lvalue = 0.0;
                  }
                  // labelLikelihoods[k - 1]->data[i] = labelLikelihoods[k - 1]->data[i] 
                  //                     + (atlasImages[j]->data[i]) * lvalue; 
                  tmpval = labelLikelihoods[k - 1].getVoxelValue(x,y,z) + atlasImages[j].getVoxelValue(x,y,z) * lvalue; 
                  labelLikelihoods[k - 1].putVoxelValue(x,y,z,tmpval);
                }
              }
            }
          }
        }
      }
    }
  }
  return(0);
}
 
 
// minimizes the posterior by the ICM algorithm. returns the number of iterations required
// labels that are forbidden for each region should receive the probability 0 in the mixture. 
// no explicit checking for forbidden labels is done


int computeMRF(LabelImage& labels,MixtureSpec& mixture,
    const FloatImage& mask,const std::vector<FloatImage>& labelLikelihoods, 
    const std::vector<FloatImage>& atlasImages, float beta1, float beta2,int maxIterations, bool verbose)
{
  int x,y,z;
  int i,j,k,r;
  char l;
  int dimx,dimy,dimz;
  int numberOfLabels,numberOfRegions;
  
 
  bool changed = true;
  int iteration = 0;
  float distLookup[27];
  float sliceWidth[3], sliceWidthMin;
  char newLabel;
  FloatImage tmpimg;

  //  bool testverbose;
  // int count = 0;

  if(verbose) cout << "Initializing... " << endl;
  // save few commonly needed values  
  dimx = mask.header.x_dim;
  dimy = mask.header.y_dim;  
  dimz = mask.header.z_dim;
  numberOfLabels  = mixture.patlas->numberOfLabels;
  numberOfRegions =  mixture.patlas->n;
  sliceWidthMin = MIN(fabs(labels.header.step[0]),fabs(labels.header.step[1]));
  sliceWidthMin = MIN(sliceWidthMin,fabs(labels.header.step[2]));

  sliceWidth[0] = (labels.header.step[0]) / sliceWidthMin;
  sliceWidth[1] = (labels.header.step[1]) / sliceWidthMin;
  sliceWidth[2] = (labels.header.step[2]) / sliceWidthMin;
  cout << sliceWidth[0] << " " << sliceWidth[1] << " " << sliceWidth[2] << " " << endl;
  // allocate
  std::vector<float> voxelProb(numberOfLabels - 1);
  std::vector<float> gibbsProb(numberOfLabels - 1);  
  std::vector<float> posteriorProb(numberOfLabels - 1);

  // fill the look up table
  for(i = (-1);i < 2;i++) {
    for(j = (-1);j < 2;j++) {
      for(k = (-1);k < 2;k++) {
        distLookup[ (i + 1) * 9 + (j + 1) * 3 + k + 1 ] = sqrt(pow(sliceWidth[0] * abs(i),2) +
                                                    pow(sliceWidth[1] * abs(j),2) +
                                                    pow(sliceWidth[2] * abs(k),2));
      }
    }
  }
  // start by initializing the ICM
  for(x = 0; x < dimx ; x++) {
    for(y = 0; y < dimy; y++) {
      for(z = 0;z < dimz; z++) {
	//         if((x == 160) && (y == 160) && (z == 190))  testverbose = true;
        // else testverbose = false;  
        if(mask.getVoxelValue(x,y,z) > 0.5) {
          collectValuesFromImagePP(labelLikelihoods,&voxelProb[0],x,y,z,numberOfLabels - 1);
          for(l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = 0.0;
          }
          for(r = 0;r < numberOfRegions;r++) {
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = posteriorProb[l]  +  
                  atlasImages[r].getVoxelValue(x,y,z) * (exp ( beta1 * log (mixture.getProb(r,(l + 1)) + 0.0001)));
            }
          }
          for(l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = posteriorProb[l] * voxelProb[l];
          }
          labels.putLabelValue(x,y,z,maxArg(&posteriorProb[0],numberOfLabels - 1) + 1);
          
	  //  count++;
	  //  if( testverbose ) cout << "test " <<  (int) (maxArg(voxelProb,numberOfLabels - 1) + 1) << endl;
        }
        else {
          labels.putLabelValue(x,y,z,0);
          // if( testverbose ) cout << "test2 " << endl; 
        }
      }
    }
  }
  

  // start the ICM
  while(changed && (iteration < maxIterations)) {
    changed = false;
    iteration++;
    if( verbose ) cout << "iteration " << iteration << endl;
    for(x = 0; x < dimx ; x++) {
      for(y = 0; y < dimy; y++) {
        for(z = 0;z < dimz; z++) {
          if(mask.getVoxelValue(x,y,z) > 0.5) {
            
            //  compute the second order term in the prior, non-normalized
            // this is the same for every region.
            for(l = 0;l < (numberOfLabels - 1);l++) {
              gibbsProb[l] = secondOrderGibbs(l + 1,labels,&mixture.patlas->mrfConstants[l + 1][0],
                                              x,y,z,distLookup,beta2);
	      //  if(testverbose) cout << " " << gibbsProb[l];
            }
	    //            if(testverbose) cout << endl;
            // then start to calculate the posterior probabilities region by region
            // initialize label probabilities to 0
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = 0.0;
            }
            // then compute region wise prior
            for(r = 0;r < numberOfRegions;r++) {
              for(l = 0;l < (numberOfLabels - 1);l++) {
                voxelProb[l] = ( gibbsProb[l] ) * (exp ( beta1 * log (mixture.getProb(r,(l + 1)) + 0.0001)));
              }
              normalize(&voxelProb[0],numberOfLabels - 1);
              for(l = 0;l < (numberOfLabels - 1);l++) {
                posteriorProb[l] = posteriorProb[l] 
                                 + atlasImages[r].getVoxelValue(x,y,z) * voxelProb[l];
              }    
            }
	    //  if(testverbose) {
            //  for(l = 0;l < (numberOfLabels - 1);l++) {
	    //   cout << " " <<  posteriorProb[l];
            //  }
            //  cout << endl;
            // } 
            // now the prior probability is computed and it remains to 
            // multiply it with the likelihood term to get the posterior
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] =  posteriorProb[l] * labelLikelihoods[l].getVoxelValue(x,y,z);
            }
            // then just find the minimum posterior and update the label
            newLabel = maxArg(&posteriorProb[0],(numberOfLabels - 1)) + 1;
            if(newLabel != labels.getLabelValue(x,y,z)) {
              changed = true;
              labels.putLabelValue(x,y,z,newLabel);  
            }
            // if(em) {
	    // cout << "Re-estimating parameters" << endl;
            //  if(mixture.patlas->labelTypes[newLabel].pureLabel) {
                

          } // endif
        }   // end for z
      }     // end for y
    }       // end for x
  }         // endwhile

  return(iteration);
}

int computeGibbs(LabelImage& labels,MixtureSpec& mixture, const FloatImage& mask,const std::vector<FloatImage>& labelLikelihoods, 
    const std::vector<FloatImage>& atlasImages, float beta1,float beta2,int maxIterations, bool verbose )
{
  int x,y,z;
  int i,j,k,r;
  char l;
  int dimx,dimy,dimz;
  int numberOfLabels,numberOfRegions;
  bool changed = true;
  int iteration = 0;
  float distLookup[27];
  float sliceWidth[3], sliceWidthMin;
  char newLabel;
  FloatImage tmpimg;

  bool testverbose;

  // int count = 0;
  if(verbose) cout << "Initializing... " << endl;
  
  // save few commonly needed values  
  dimx = mask.header.x_dim;
  dimy = mask.header.y_dim;  
  dimz = mask.header.z_dim;
  numberOfLabels = mixture.patlas->numberOfLabels;
  numberOfRegions =  mixture.patlas->n;
  sliceWidthMin = MIN(fabs(labels.header.step[0]),fabs(labels.header.step[1]));
  sliceWidthMin = MIN(sliceWidthMin,fabs(labels.header.step[2]));

  sliceWidth[0] = (labels.header.step[0]) / sliceWidthMin;
  sliceWidth[1] = (labels.header.step[1]) / sliceWidthMin;
  sliceWidth[2] = (labels.header.step[2]) / sliceWidthMin;
  cout << sliceWidth[0] << " " << sliceWidth[1] << " " << sliceWidth[2] << " " << endl;
  // allocate
  
  std::vector<float> voxelProb(numberOfLabels - 1);
  std::vector<float> gibbsProb(numberOfLabels - 1);  
  std::vector<float> posteriorProb(numberOfLabels - 1);

  // fill the look up table
  for(i = (-1);i < 2;i++) {
    for(j = (-1);j < 2;j++) {
      for(k = (-1);k < 2;k++) {
        distLookup[ (i + 1) * 9 + (j + 1) * 3 + k + 1 ] = sqrt(pow(sliceWidth[0] * abs(i),2) +
                                                    pow(sliceWidth[1] * abs(j),2) +
                                                    pow(sliceWidth[2] * abs(k),2));
      }
    }
  }
  // start by initializing the ICM
  for(x = 0; x < dimx ; x++) {
    for(y = 0; y < dimy; y++) {
      for(z = 0;z < dimz; z++) {
        if((x == 190) && (y == 32) && (z == 105))  testverbose = true;
        else testverbose = false;  
        if(mask.getVoxelValue(x,y,z) > 0.5) {
          collectValuesFromImagePP(labelLikelihoods,&voxelProb[0],x,y,z,numberOfLabels - 1);
          for(l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = 0.0;
          }
          for(r = 0;r < numberOfRegions;r++) {
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = posteriorProb[l]  +  atlasImages[r].getVoxelValue(x,y,z) * 
                                 (exp ( beta1 * log (mixture.getProb(r,(l + 1)) + 0.0001)));
            }
          }
          if(testverbose) {
            for(l = 0;l < (numberOfLabels - 1);l++) {
              cout << " " << voxelProb[l];
            } 
            cout << endl;
             for(l = 0;l < (numberOfLabels - 1);l++) {
              cout << " " << posteriorProb[l];
            }    
            cout << endl;
             for(r = 0;r < numberOfRegions;r++) {
              cout << " " << atlasImages[r].getVoxelValue(x,y,z);
            } 
	     cout << endl;
          }
          for(l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = posteriorProb[l] * voxelProb[l]; // MAP init
	    //  posteriorProb[l] = voxelProb[l]; // mlinit 
          }
          labels.putLabelValue(x,y,z,maxArg(&posteriorProb[0],numberOfLabels - 1) + 1);
	  //  count++;
	  //  if( testverbose ) cout << "test " <<  (int) (maxArg(voxelProb,numberOfLabels - 1) + 1) << endl;
        }
        else {
          labels.putLabelValue(x,y,z,0);
          // if( testverbose ) cout << "test2 " << endl; 
        }
      }
    }
  }
  
  // writeLabelImage("tmp",labels,true);
  // start the ICM
  while(changed && (iteration < maxIterations)) {
    changed = false;
    iteration++;
    if( verbose ) cout << "iteration " << iteration << endl;
    for(x = 0; x < dimx ; x++) {
      for(y = 0; y < dimy; y++) {
        for(z = 0;z < dimz; z++) {
          if(mask.getVoxelValue(x,y,z) > 0.5) {
            
            //  compute the second order term in the prior, non-normalized
            // this is the same for every region.
            for(l = 0;l < (numberOfLabels - 1);l++) {
              gibbsProb[l] = secondOrderGibbs(l + 1,labels,&mixture.patlas->mrfConstants[l + 1][0],
                                              x,y,z,distLookup,beta2);
	      //  if(testverbose) cout << " " << gibbsProb[l];
            }
	    //            if(testverbose) cout << endl;
            // then start to calculate the posterior probabilities region by region
            // initialize label probabilities to 0
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = 0.0;
            }
            // then compute region wise prior
            for(r = 0;r < numberOfRegions;r++) {
              for(l = 0;l < (numberOfLabels - 1);l++) {
                voxelProb[l] = ( gibbsProb[l] ) * (exp ( beta1 * log (mixture.getProb(r,(l + 1)) + 0.0001)));
              }
	      //  normalize(voxelProb,numberOfLabels - 1); this is the difference between Markov and Gibbs
              for(l = 0;l < (numberOfLabels - 1);l++) {
                posteriorProb[l] = posteriorProb[l] 
                                 + atlasImages[r].getVoxelValue(x,y,z) * voxelProb[l];
              }    
            }
            
	    //  if(testverbose) {
            //  for(l = 0;l < (numberOfLabels - 1);l++) {
	    //   cout << " " <<  posteriorProb[l];
            //  }
            //  cout << endl;
            // } 
            // now the prior probability is computed and it remains to 
            // multiply it with the likelihood term to get the posterior
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] =  posteriorProb[l] * labelLikelihoods[l].getVoxelValue(x,y,z);
            }
            // then just find the minimum posterior and update the label
            newLabel = maxArg(&posteriorProb[0],(numberOfLabels - 1)) + 1;
            if(newLabel != labels.getLabelValue(x,y,z)) {
              changed = true;
              labels.putLabelValue(x,y,z,newLabel);  
            }
          } // endif
        }   // end for z
      }     // end for y
    }       // end for x
  }         // endwhile

  return(iteration);
}

int convertPVElabels(LabelImage& crispLabels, const LabelImage& pveLabels, const FloatImage& img, 
    const std::vector<FloatImage>& atlasImages, MixtureSpec& mixture)
{
  int n,r;
  int i,j,k;
  char label,pureLabel1,pureLabel2;;
  float likelihoodValue[INTERVALS + 1];
  float rprob,mean1,mean2,var1,var2;;
  float maxindex,maxvalue;
  float t;

  for(i = 0;i < pveLabels.header.x_dim;i++) {
    for(j = 0;j < pveLabels.header.y_dim;j++) {
      for(k = 0;k < pveLabels.header.z_dim;k++) {
        label = pveLabels.getLabelValue(i,j,k);
        if(mixture.patlas->labelTypes[label].pureLabel) {
          crispLabels.putLabelValue(i,j,k,label);
        }
        else {
          pureLabel1 = mixture.patlas->labelTypes[label].mixed[0]; 
          pureLabel2 = mixture.patlas->labelTypes[label].mixed[1];
          if(pureLabel1 == 0) {
             crispLabels.putLabelValue(i,j,k,pureLabel2);
          }
          else if(pureLabel2 == 0) {
            crispLabels.putLabelValue(i,j,k,pureLabel1);
          }
          else {
            for(n = 0;n < (INTERVALS + 1);n++) {
              likelihoodValue[n] = 0.0;
            }
            for(r = 0;r < mixture.patlas->n;r++) {
              rprob = atlasImages[r].getVoxelValue(i,j,k);
              if(rprob > 0.001) {
                mean1 = mixture.getMu(r,pureLabel1);
                mean2 = mixture.getMu(r,pureLabel2);          
                var1 = mixture.getSigma2(r,pureLabel1);
                var2 = mixture.getSigma2(r,pureLabel2);
                for(n = 0;n < (INTERVALS + 1);n++) {
                  t = n * (float) 1/ INTERVALS;
                  likelihoodValue[n] = likelihoodValue[n] + rprob * 
		                       computeGaussianLikelihood(img.getVoxelValue(i,j,k),
                                                               t * mean1 + ( 1 - t) * mean2,
                                                               t * t * var1 + (1 - t) * ( 1-  t) *var2);
                }
              }
            }
            maxvalue = likelihoodValue[0];
            maxindex = 0;
            for(n = 0;n < (INTERVALS + 1) ;n++) {
              if(likelihoodValue[n] > maxvalue) {
                maxvalue = likelihoodValue[n];
                maxindex = n;
              }
            } 
            if(maxindex > (INTERVALS/2)) {
              crispLabels.putLabelValue(i,j,k,pureLabel1);   
            }   
            else {
              crispLabels.putLabelValue(i,j,k,pureLabel2);
            }
          } // end else
        } // end else
      }   // end for k
    }     // end for j
  }       // end for i
  return(0);
}
