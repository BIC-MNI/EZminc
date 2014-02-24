/*
 *  dwt.h
 *  Wavelet
 *
 *  Created by Robert Brown on 26/05/09.
 *  Copyright 2009 University of Calgary. All rights reserved.
 *
 */

#include <minc_io_simple_volume.h>

namespace minc
{

	struct wavefilt { 
		int ncof,ioff,joff; 
		float *cc,*cr; 
	} ; 


	void daub4(float a[], int n, int stride, int isign, wavefilt wfilt);
	void pwt(float a[],  int n, int stride, int isign, wavefilt wfilt);

	void dwt(float a[], int n, int top, int stride, int isign, void (*wtstep)(float [], int, int, int, wavefilt));

	void dwt2DNonStandard(float a[], int N[3], int n[2], int top, int stride, int isign, void (*wtstep)(float [], int, int, int, wavefilt));

	void dwt3DNonStandard(float* a, int* N, int* n, int top, unsigned int stride, int isign, void (*wtstep)(float [], int, int, int, wavefilt));

	void volume_dwt(simple_volume<float>& vol,bool forward,void (*wtstep)(float [], int, int, int, wavefilt)=pwt); 

};
