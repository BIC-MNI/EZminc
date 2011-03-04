#ifndef _MINC_WRAPPERS_H_
#define _MINC_WRAPPERS_H_

#include <complex>
#include <vector>
#include <algorithm>
#include <itkArray.h>
#include <iostream>

#include <minc_io_exceptions.h>
#include <minc_io_fixed_vector.h>
#include <minc_helpers.h>

namespace minc
{
  //! minc4itk compatibility function
  template <class T> void  load_minc(const char *file,T &img)
  {
    img=load_minc<typename T::ObjectType>(file);
  }
  /*
  //! minc4itk compatibility function
  template <class T> void  save_minc(const char *file,const T &img)
  {
    save_minc<typename T::ObjectType>(file,img);
  }*/
  
  template<class T> void setup_itk_image(const minc_1_base& minc_rw, T& img)
  {
    itk::Vector< unsigned int,3> dims;
    itk::Vector< double,3> spacing;
    itk::Vector< double,3> origin;
    itk::Vector< double,3> start;
    itk::Matrix< double, 3, 3> dir_cos;
    dir_cos.SetIdentity();
    //std::cout<<"setup_itk_image"<<std::endl;
    for(int i=0;i<3;i++)
    {
      dims[i]=minc_rw.ndim(i+1);
      spacing[i]=minc_rw.nspacing(i+1);
      start[i]=minc_rw.nstart(i+1);
      if(minc_rw.have_dir_cos(i+1))
      {
        for(int j=0;j<3;j++)
          dir_cos[j][i]=minc_rw.ndir_cos(i+1,j); //TODO: check transpose?
      }
      //std::cout<<start[i]<<"\t";
    }
    //std::cout<<std::endl;
    origin=dir_cos*start;
    allocate_image3d(img,dims,spacing,origin);
    img->SetDirection(dir_cos);
  }
  
  template<class T> void imitate_minc (const char *path, T& img)
  {
    minc_1_reader rdr;
    rdr.open(path,true,true);
    setup_itk_image<T>(rdr,img);
  }
  
  
};//minc
#endif //_MINC_WRAPPERS_H_
