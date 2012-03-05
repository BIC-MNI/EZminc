#ifndef __MINC_GENERIC_IMAGE_H__
#define __MINC_GENERIC_IMAGE_H__

#include "generic_image.h"

int readImage(const char* filename,  FloatImage& img);
int writeImage(const char* filename,const FloatImage& img);
int readImage(const char* filename,LabelImage& img);
int writeImage(const char* filname,const LabelImage& img);


#endif //__MINC_GENERIC_IMAGE_H__
