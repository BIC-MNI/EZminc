/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2006 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include "dwt_utils.h"
#include "dwt.h"
#include <assert.h>

namespace minc
{
  idx3d find_nearest_square_pow2(const idx3d& i)
  {
    unsigned int m=i.max(),_m=i.max();
    unsigned int p=0;
  
    while(m!=0) { m>>=1;p++;}
  
    if(_m==(1<<(p-1))) p--; //special case, when we have exatly power of 2
  
    return IDX(1<<p,1<<p,1<<p);
  }

  void dwt_forward(const simple_volume<float> &src,std::vector<simple_volume<float> > &dst)
  {
    idx3d input_size;
    idx3d padded_size;
    idx3d output_size;
    idx3d pad;
    
    input_size=src.size();
    
    output_size=(input_size+IDX(1,1,1))/2;
    
    padded_size=find_nearest_square_pow2(input_size);
    pad=(padded_size-input_size)/2;
    
    simple_volume<float> input_vol;
    input_vol.resize(padded_size);
    
    minc::pad_volume<float>(src,input_vol,0); //pad volume with zeros
    
    volume_dwt(input_vol,1);
    
    pad/=2;
    dst.resize(8);
    for(int j=0;j<8;j++)
    {
      dst[j].resize(output_size);
      extract_quadrant<float>(input_vol,dst[j],j,pad);
    }
  }
  
  void dwt_backward(const std::vector<simple_volume<float> > &src,simple_volume<float> &dst)
  {
    assert(src.size()==8);
    
    idx3d input_size;
    idx3d padded_size;
    idx3d output_size;
    idx3d pad;
    
    input_size=src[0].size();
    output_size=input_size*2;
    
    padded_size=find_nearest_square_pow2(input_size);
    
    pad=(padded_size-input_size)/2;
    
    simple_volume<float> input_vol(padded_size[0]*2,padded_size[1]*2,padded_size[2]*2);
    
    for(int j=0;j<8;j++)
      insert_quadrant<float>(input_vol,src[j],j,pad);
      
    volume_dwt(input_vol,0);
    
    dst.resize(output_size);
    
    minc::pad_volume<float>(input_vol,dst,0); 
  }


};
