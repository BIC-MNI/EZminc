/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2009 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <string>
#include <minc_io_simple_volume.h>

namespace minc
{
  
  //! vector of strings 
  typedef std::vector<std::string> strings;
  
  //! vector of doubles 
  typedef std::vector<double> doubles;
  
  //! string table (vector of vectors of strings)
  typedef std::vector<strings> string_table;
  
  //! doubles table (vector of vectors of doubles)
  typedef std::vector<doubles> doubles_table;
  
  //! vector of minc volumes
  typedef std::vector<minc_float_volume > volumes;
  //! vector of minc grids
  typedef std::vector<minc_grid_volume > grids;
  
  //! convert list of comma separated values into vector
  void string_to_vector(const std::string& s,std::vector<double>& l,const char *ws=" ,");
  
  //! read a volume from a csv file (without a header), lines starting with # are commments
  void read_table(const char *fn,string_table& tbl);
  
  //! calculate column-wise averages
  void average_tables(string_table& tbl,volumes &v);
  //! calculate column-wise averages
  void average_tables(string_table& tbl,grids &v);
  
  //! make sure that table is square 
  bool table_verify(const string_table& tbl);
  
  //! make sure that table is square 
  bool table_verify(const doubles_table& tbl);
  
  //! read a volume from a csv file (without a header), lines starting with # a commments, entries with colums > limit are trated as relative file names, so the directory of the table is added in front
  void read_table_n(const char *fn,string_table& tbl,int limit=1);
  
  //! read a volume from a csv file (without a header), lines starting with # a commments, entries with colums < limit are trated as relative file names, so the directory of the table is added in front
  void read_table_f(const char *fn,string_table& tbl,int limit=1);

  //! read a doubles table, lines starting with # are comments
  void read_doubles_table(const char *fn,doubles_table& tbl);
  
  //! calculate volume variance
  double variance(const minc_float_volume& v,const minc_float_volume& mean);
  
  //! calculate sum of all voxels in volume
  double sum(const minc_float_volume& v);
  //! calculate sum of all voxels in volume
  double sum(const minc_grid_volume& v);
  //! calculate sum of all voxels in volume
  double sum(const minc_byte_volume& v);
  
  //! calculate sum of all voxels^2 in volume, withing the mask
  double sum2(const minc_float_volume& v,const minc_byte_volume& mask);
  //! calculate sum of all voxels^2 in volume
  double sum2(const minc_float_volume& v);
  //! calculate sum of all ||vectors||^2 in volume, withing the mask
  double sum2(const minc_grid_volume& v,const minc_byte_volume& mask);
  
  //! calculate product of two images 
  void mul(minc_float_volume &v2,const minc_float_volume& v1);
  
  //! calculate product of two images within mask
  void masked_mul(minc_float_volume &v2,const minc_float_volume& v1,const minc_byte_volume &mask);
  
  //! calculate product of two grid volumes within mask
  void masked_mul(minc_grid_volume &v2,const minc_grid_volume& v1,const minc_byte_volume &mask);
  
  //! load a minc volume
  void load_volume(const char *fn,minc::minc_grid_volume& vol);
  //! load a minc volume
  void load_volume(const char *fn,minc::minc_byte_volume& vol);
  //! load a minc grid
  void load_volume(const char *fn,minc::minc_float_volume& vol);
  
  //! load all volumes , starting from skip
  void load_volumes(strings& s,volumes& vols,int skip=1,bool verbose=false,bool ignore_missing=false);
  //! load all grids , starting from skip
  void load_volumes(strings& s,grids& vols,int skip=1,bool verbose=false);
  //! check if volumes are same dimensions
  void check_volumes(volumes& vols,bool ignore_missing=false);
  //! save a minc volume
  void save_volume(minc::minc_grid_volume& vol,const char *fn,const char *like,const std::string & append_history="");
  //! save a minc volume
  void save_volume(minc::minc_byte_volume& vol,const char *fn,const char *like,const std::string & append_history="");
  //! save a minc grid
  void save_volume(minc::minc_float_volume& vol,const char *fn,const char *like,const std::string & append_history="");
  //! save a minc volumes
  void save_volumes(strings& s,volumes& vols,const char *like,const std::string & append_history="");
  
  
  //!load all basis volumes, skipping the first line (usuallt the mean)
  void load_all_volumes(string_table & tbl,std::vector<volumes>& all_volumes,int selected,bool verbose=false);
  
  //!load all basis volumes
  void load_all_volumes_0(string_table & tbl,std::vector<volumes>& all_volumes,int selected,bool verbose=false);
  
  //!load all basis grids, skipping the first line (usuallt the mean)
  void load_all_volumes(string_table & tbl,std::vector<grids>& all_volumes,int selected);

  //! v1-=v2
  void substract(volumes& v1,const volumes& v2);
  //! v1-=v2
  void substract(grids& v1,const grids& v2);
  //! analog to dirname command
  std::string  _dirname(const char *file);
  //! analog to basename command
  std::string _basename(const char *file);
  
  //! calculate voxel volume 
  double calc_dv(const char* minc_file);
  
  //! calculate RMS difference between volumes
  double volume_rms_diff(const minc::minc_float_volume& in1,const minc::minc_float_volume& in2);

  
  
};
#endif //__UTILS_H__
