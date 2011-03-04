/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: interface to 3D wavelet transform
@COPYRIGHT  :
              Copyright 2009 Robert Brown,Vladimir Fonov, 
              McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include <minc_io_simple_volume.h>
#include <minc_io_fixed_vector.h>
#include <vector>

namespace minc
{
  typedef fixed_vec<3,int> idx3d;
  
  template<class T> void extract_quadrant(const simple_volume<float> &src,simple_volume<float>& dst,int quadrant,const idx3d& pad=IDX<int>(0,0,0))
  {
//make sure out is one half
    for(int i=0;i<3;i++)
      if((src.dim(i)/2)<((dst.dim(i)+pad[i])))
        REPORT_ERROR("Unexpected dimension size");
      
    int qx=( quadrant&1    )*src.dim(0)/2+pad[0];
    int qy=((quadrant&2)>>1)*src.dim(1)/2+pad[1];
    int qz=((quadrant&4)>>2)*src.dim(2)/2+pad[2];
  
    for(int z=0;z<dst.dim(2);z++)
      for(int y=0;y<dst.dim(1);y++)
        for(int x=0;x<dst.dim(0);x++)
    {
      dst.set(x,y,z,src.get(x+qx,y+qy,z+qz));
    }
  }

  template<class T> void insert_quadrant(simple_volume<float> &dst,const simple_volume<float>& src,int quadrant,const idx3d& pad=IDX<int>(0,0,0))
  {
    for(int i=0;i<3;i++)
      if((dst.dim(i)/2)<((src.dim(i)+pad[i])))
        REPORT_ERROR("Unexpected dimension size");
      
    int qx=(quadrant&1)     *dst.dim(0)/2+pad[0];
    int qy=((quadrant&2)>>1)*dst.dim(1)/2+pad[1];
    int qz=((quadrant&4)>>2)*dst.dim(2)/2+pad[2];
  
    for(int z=0;z<src.dim(2);z++)
      for(int y=0;y<src.dim(1);y++)
        for(int x=0;x<src.dim(0);x++)
    {
      dst.set(x+qx,y+qy,z+qz,src.get(x,y,z));
    }
  }
  
  idx3d find_nearest_square_pow2(const idx3d&i);
  
  
  void dwt_forward(const simple_volume<float> &src,std::vector<simple_volume<float> > &dst);
  void dwt_backward(const std::vector<simple_volume<float> > &src,simple_volume<float> &dst);
  
};
