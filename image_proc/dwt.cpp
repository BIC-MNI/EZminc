/*
 *  dwt.c
 *  Wavelet
 *
 *  Created by Robert Brown on 26/05/09.
 *  Copyright 2009 University of Calgary. All rights reserved.
 *
 *  Modified by Vladimir Fonov on 25/09/09
 */

#include <stdlib.h>
#include <math.h>
#include "dwt.h"


#define C0 0.4829629131445341 
#define C1 0.8365163037378079 
#define C2 0.2241438680420134 
#define C3 -0.1294095225512604 

namespace minc
{

	// One level wavelet decomposition using the Daubechies D4 wavelet
	void daub4(float a[], int n, int stride, int isign, wavefilt wfilt) {
		float *wksp;
		int nh, nh1, i, j;
		
		if (n < 4) return;
		wksp = (float *)malloc(sizeof(float) * n);
		nh1 = (nh=n >> 1) + 1;
		if (isign >= 0) {
			for (i=0, j=0; j<n-3; j+=2, i++) {
				wksp[i] = C0*a[j*stride] + C1*a[(j+1)*stride] + C2*a[(j+2)*stride] + C3*a[(j+3)*stride];
				wksp[i+nh] = C3*a[j*stride] - C2*a[(j+1)*stride] + C1*a[(j+2)*stride] - C0*a[(j+3)*stride];
			}
			wksp[i] = C0*a[(n-2)*stride] + C1*a[(n-1)*stride] + C2*a[0] + C3*a[1*stride];
			wksp[i+nh] = C3*a[(n-2)*stride] - C2*a[(n-1)*stride] + C1*a[0] - C0*a[1*stride];
		} else {
			wksp[0] = C2*a[(nh-1)*stride] + C1*a[(n-1)*stride] + C0*a[0] + C3*a[(nh1-1)*stride];
			wksp[1] = C3*a[(nh-1)*stride] - C0*a[(n-1)*stride] + C1*a[0] - C2*a[(nh1-1)*stride];
			for (i=0, j=2; i<nh-1; i++) {
				wksp[j++] = C2*a[i*stride] + C1*a[(i+nh)*stride] + C0*a[(i+1)*stride] + C3*a[(i+nh1)*stride];
				wksp[j++] = C3*a[i*stride] - C0*a[(i+nh)*stride] + C1*a[(i+1)*stride] - C2*a[(i+nh1)*stride];
			} 
		}
		for (i=0; i<n; i++)
			a[i*stride]=wksp[i];
		free(wksp);
	}

	void pwt(float a[],  int n, int stride, int isign, wavefilt wfilt) { 
		float ai, ai1, *wksp; 
		int i, ii, j, jf, jr, k, n1, ni, nj, nh, nmod; 
		if (n < 4) return; 
		wksp = (float *) malloc(sizeof(float)*n); 
		nmod = wfilt.ncof*n;
		n1 = n-1;
		nh = n >> 1;
		for (j=0; j<n; j++)
			wksp[j] = 0.0; 
		if (isign >= 0) {
			for (ii=0,i=0;i<n;i+=2,ii++) { 
				ni = i+nmod+wfilt.ioff;
				nj = i+nmod+wfilt.joff; 
				for (k=0; k<wfilt.ncof; k++) { 
					jf = (n1&(ni+k+2))-1;
					jr = (n1&(nj+k+2))-1; 
					wksp[ii] += wfilt.cc[k]*a[(jf+1)*stride]; 
					wksp[ii+nh] += wfilt.cr[k]*a[(jr+1)*stride]; 
				} 
			} 
		} else {
			for (ii=0,i=0;i<n;i+=2,ii++) { 
				ai = a[ii*stride]; 
				ai1 = a[(ii+nh)*stride]; 
				ni = i+nmod+wfilt.ioff;
				nj = i+nmod+wfilt.joff; 
				for (k=0;k<wfilt.ncof;k++) { 
					jf = (n1&(ni+k+2)); 
					jr = (n1&(nj+k+2)); 
					wksp[jf] += wfilt.cc[k]*ai; 
					wksp[jr] += wfilt.cr[k]*ai1; 
				} 
			} 
		} 
		for (j=0;j<n;j++)
			a[j*stride]=wksp[j];
		free(wksp);
	}

	void dwt(float a[], int n, int top, int stride, int isign,
			void (*wtstep)(float [], int, int, int, wavefilt)) {
		int nn;
		static float cc[4] = {C0,C1,C2,C3};
		static float cr[4] = {C3,-C2,C1,-C0};
		wavefilt wfilt;
		wfilt.cc = cc;
		wfilt.cr = cr;
		wfilt.ncof = 4;
		wfilt.ioff = wfilt.joff = -2;

		if (n < 4) return;
		if (top) {
			if (isign >= 0) {
				for (nn=n; nn>=4; nn>>=1)
					(*wtstep)(a,nn,stride,isign,wfilt);
			} else {
				for (nn=4; nn<=n; nn<<=1)
					(*wtstep)(a,nn,stride,isign,wfilt);
			}
		} else {
			if (isign >= 0) {
				(*wtstep)(a,n,stride,isign,wfilt);
			} else {
				(*wtstep)(a,n,stride,isign,wfilt);
			}
		}
	}


	void dwt2DNonStandard(float a[], int N[2], int n[2], int top, int stride, int isign,
						void (*wtstep)(float [], int, int, int, wavefilt)) {
		int startStride;
		unsigned int lstride;
		int done = 0;
		int iterated[2];
		int norig[2];
		
		norig[0] = n[0];
		norig[1] = n[1];

		if (isign <= 0 && top) {
			n[0] = 4;
			n[1] = 4;
		}

		while (!done) {
			iterated[1] = 0;
			if (n[1] >= n[0]) {
				// All the rows
				startStride = N[1]*stride;
				lstride = stride;
				if (n[1] <= N[1]) {
					for (int i=0; i < fmin(n[0],N[0]); i++) {
						dwt(a+i*startStride, n[1], 0, lstride, isign, wtstep);
					}
					iterated[1] = 1;
				}
			}
			iterated[0] = 0;		
			if (n[0] >= n[1]) {
				// All the cols
				lstride = N[1]*stride;
				startStride = stride;
				if (n[0] <= N[0]) {
					for (int i=0; i < fmin(n[1],N[1]); i++) {		
						dwt(a+i*startStride, n[0], 0, lstride, isign, wtstep);
					}
					iterated[0] = 1;
				}
			}

			if (isign > 0) {
				if (iterated[0])
					n[0] = n[0] >> 1;
				if (iterated[1])
					n[1] = n[1] >> 1;
				if ((n[1] < 4 && n[0] < 4) || !top)
					done = 1;
			} else {
				if (iterated[0])
					n[0] = n[0] << 1;
				if (iterated[1])
					n[1] = n[1] << 1;
				if ((n[1] > N[1] && n[0] > N[0]) || !top)
					done = 1;
			}
		}
		n[0] = norig[0];
		n[1] = norig[1];
	}

	void dwt3DNonStandard(
	float* a, int* N, int* n, int top, unsigned int stride, int isign,
							void (*wtstep)(float [], int, int, int, wavefilt)) {
		int startStride[2];
		unsigned int lstride;
		int done = 0;
		int iterated[3];
		int doneDims[3];

		for (int i=0; i<3; i++) {
			n[i] = N[i];
			doneDims[i] = 0;
		}
	/*
		if (isign <= 0) {
			n[0] = N[0]/2;
			n[1] = N[1]/2;
			n[2] = N[2]/2;
		}	*/
		
		while (!done) {
		
			for (int i=0; i<3; i++) {
				if (n[i] > N[i]) {
					doneDims[i] = 1;
	//				n[i] = N[i];
				}
			}
			// In plane
			iterated[1] = 0;
			iterated[2] = 0;
			if (n[1] >= n[2] && n[1] >= n[0] && isign > 0)
				iterated[1] = 1;
			if (n[1] <= N[1] && isign <= 0)
				iterated[1] = 1;			
			if (n[2] >= n[1] && n[2] >= n[0])
				iterated[2] = 1;
			if (n[2] <= N[2] && isign <= 0)
				iterated[2] = 1;			
			startStride[0] = N[2]*N[1]*stride;
			lstride = stride;
			//fprintf(stderr,"Transforming plane with n=[%d,%d]\n",n[1],n[2]);
			for (int i=0; i < fmin(n[0],N[0]); i++) {
				dwt2DNonStandard(a+i*startStride[0], N+1, n+1, 0, lstride,isign,wtstep);
			}

			// through plane
			iterated[0] = 0;
			if ((n[0] >= n[1]) && (n[0] >= n[2])) {		
				lstride = N[2]*N[1]*stride;
				startStride[0] = stride;
				startStride[1] = stride * N[2];
				if (!doneDims[0]) {			
					//fprintf(stderr,"Transforming through planes for %d x %d x %d\n",n[0],n[1],n[2]);
					for (int h=0; h < fmin(n[2],N[2]); h++) {
						for (int v=0; v < fmin(n[1],N[1]); v++ ) {
	//						if (i % 4 == 0) fprintf(stderr,"Doing dwt for n=%d starting at %d with stride %d\n",n[0],i*startStride,lstride);
							dwt(a+h*startStride[0]+v*startStride[1], n[0], 0, lstride, isign, wtstep);
						}
					}
					iterated[0] = 1;
				}
			}
			
			
			if (isign > 0) {
				int steps = 4;
				if (iterated[0])
					n[0] = n[0] >> 1;
				if (iterated[1])
					n[1] = n[1] >> 1;
				if (iterated[2])
					n[2] = n[2] >> 1;
				if ((n[2] < steps && n[1] < steps && n[0] < steps) || !top)
					done = 1;
			} else {
				if (iterated[0])
					n[0] = n[0] << 1;
				if (iterated[1])
					n[1] = n[1] << 1;
				if (iterated[2])
					n[2] = n[2] << 1;
				if (((n[2] > N[2] || !iterated[2]) && (n[1] > N[1] || !iterated[1]) && (n[2] > N[2] || !iterated[1])) || !top)
					done = 1;
			}
			
		}
	}
	
  void volume_dwt(simple_volume<float>& vol,bool forward,void (*wtstep)(float [], int, int, int, wavefilt))
  {
    int N[3]={vol.dim(0),vol.dim(1),vol.dim(2)};
    int n[3]={vol.dim(0),vol.dim(1),vol.dim(2)};

    dwt3DNonStandard(vol.c_buf(),N,n,0,1,forward?1:-1,wtstep);
  }


};
