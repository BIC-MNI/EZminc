#include <minc_1_simple.h>
#include "minc_generic_image.h"
#include <iostream>

using namespace std;

using namespace minc;

void minc2header(const minc_1_base& minc,ImageHeader &hdr)
{
  if(minc.dim_no()!=3)
    REPORT_ERROR("Only 3D minc files are supported!");
  hdr.dims=3;
  hdr.x_dim=minc.ndim(1);
  hdr.y_dim=minc.ndim(2);
  hdr.z_dim=minc.ndim(3);
  for(int i=0;i<3;i++)
  {
    hdr.start[i]=minc.nstart(i+1);
    hdr.step[i]=minc.nspacing(i+1);
    for(int j=0;j<3;j++)
      hdr.dir_cos[i][j]=minc.ndir_cos(i+1,j);
  }
}

void header2minc(const ImageHeader &hdr,minc_info &minc)
{
  minc.clear();
  minc.resize(3);
  
  minc[2].length=hdr.x_dim;
  minc[2].dim=dim_info::DIM_X;
  
  minc[1].length=hdr.y_dim;
  minc[1].dim=dim_info::DIM_Y;
  
  minc[0].length=hdr.z_dim;
  minc[0].dim=dim_info::DIM_Z;
  
  for(int i=0;i<3;i++)
  {
    minc[2-i].start=hdr.start[i];
    minc[2-i].step =hdr.step[i];
    minc[2-i].have_dir_cos=true;
    for(int j=0;j<3;j++)
      minc[2-i].dir_cos[j]=hdr.dir_cos[i][j];
  }
  
}

int readImage(const char* filename,FloatImage& img)
{
  try
  {
    minc_1_reader rdr;
    rdr.open(filename);
    rdr.setup_read_float();
    unsigned long size=1;
    for(int i=0;i<rdr.dim_no();i++)
      size*=rdr.dim(i).length;
    
    minc2header(rdr,img.header);
    
    img.data.resize(size);
    load_standard_volume<float>(rdr,&img.data[0]);
    cout<<"Read float:"<<filename<<" done:"<<img.header.x_dim<<"x"<<img.header.y_dim<<"x"<<img.header.z_dim<<endl;
    
  } catch (const minc::generic_error & err) {
    std::cerr<<err.msg();
    return 1;
  }
  return(0);
}

int writeImage(const char* filename,const FloatImage& img) 
{
  try
  {
    minc_1_writer wrt;
    minc_info info;
    header2minc(img.header,info);
    wrt.open(filename,info,2,NC_FLOAT,1); //todo: do we need to save in floats?
    unsigned long size=1;
    for(int i=0;i<wrt.dim_no();i++)
      size*=wrt.dim(i).length;
    wrt.setup_write_float();
    save_standard_volume<float>(wrt,&img.data[0]);
  } catch (const minc::generic_error & err) {
    std::cerr<<err.msg();
    return 1;
  }
  return(0);
}

int readImage(const char* filename,LabelImage& img)
{
  try
  {
    minc_1_reader rdr;
    rdr.open(filename);
    rdr.setup_read_byte();
    unsigned long size=1;
    for(int i=0;i<rdr.dim_no();i++)
      size*=rdr.dim(i).length;
    
    minc2header(rdr,img.header);
    
    img.data.resize(size);
    load_standard_volume<unsigned char>(rdr,&img.data[0]);
    cout<<"Read label:"<<filename<<" done:"<<img.header.x_dim<<"x"<<img.header.y_dim<<"x"<<img.header.z_dim<<endl;
  } catch (const minc::generic_error & err) {
    std::cerr<<err.msg();
    return 1;
  }
  return(0);
}

int writeImage(const char* filename,const LabelImage& img)
{
  try
  {
    minc_1_writer wrt;
    minc_info info;
    header2minc(img.header,info);
    wrt.open(filename,info,2,NC_BYTE); //todo: do we need to save in floats?
    unsigned long size=1;
    for(int i=0;i<wrt.dim_no();i++)
      size*=wrt.dim(i).length;
    wrt.setup_write_byte();
    save_standard_volume<unsigned char>(wrt,&img.data[0]);
  } catch (const minc::generic_error & err) {
    std::cerr<<err.msg();
    return 1;
  }
  return(0);
}
