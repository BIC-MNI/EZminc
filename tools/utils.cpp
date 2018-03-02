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
#include "utils.h"

#include <fstream>

#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include <math.h>

#include <stdlib.h>
#include <libgen.h> //for dirname

namespace minc
{

  void read_table(const char *fn,string_table& tbl)
  {
    tbl.clear();
    //read input file line by line, ignoring empty lines and lines starting with #
    std::ifstream in(fn);
    if(!in.good())
      REPORT_ERROR("Can't open file for reading!");
    std::string dir=minc::_dirname(fn);
    std::cout<<"Dir="<<dir.c_str()<<std::endl;
    dir+="/";
    while(in.good() && !in.eof())
    {
      char tmp[32000];
      in.getline(tmp,sizeof(tmp));
      if(!strlen(tmp)) continue;
      if(tmp[0]=='#') continue;
      char *str,*saveptr;
      strings _s;
      
      for(str=tmp;;str=NULL)
      {
        const char *token = strtok_r(str, ",", &saveptr);
        if(!token || !strlen(token)) break;
        _s.push_back(dir+token);
      }
      if(_s.empty()) continue;
      tbl.push_back(_s);
    }
  }
  
  
  void read_table_n(const char *fn,string_table& tbl,int limit)
  {
    tbl.clear();
    //read input file line by line, ignoring empty lines and lines starting with #
    std::ifstream in(fn);
    if(!in.good())
      REPORT_ERROR("Can't open file for reading!");
    
    std::string dir=minc::_dirname(fn);
    dir+="/";
    while(in.good() && !in.eof())
    {
      char tmp[32000];
      in.getline(tmp,sizeof(tmp));
      if(!strlen(tmp)) continue;
      if(tmp[0]=='#') continue;
      char *str,*saveptr;
      strings _s;
      //std::cout<<tmp<<std::endl;
      int k=0;
      for(str=tmp;;str=NULL)
      {
        const char *token = strtok_r(str, ",", &saveptr);
        if(!token || !strlen(token)) break;
        if(k>=limit)
          _s.push_back(dir+token);
        else
          _s.push_back(token);
        k++;
      }
      if(_s.empty()) continue;
      tbl.push_back(_s);
    }
  }
  
  void read_table_f(const char *fn,string_table& tbl,int limit)
  {
    tbl.clear();
    //read input file line by line, ignoring empty lines and lines starting with #
    std::ifstream in(fn);
    if(!in.good())
      REPORT_ERROR("Can't open file for reading!");
    std::string dir=minc::_dirname(fn);
    dir+="/";
    while(in.good() && !in.eof())
    {
      char tmp[32000];
      in.getline(tmp,sizeof(tmp));
      if(!strlen(tmp)) continue;
      if(tmp[0]=='#') continue;
      char *str,*saveptr;
      strings _s;
      //std::cout<<tmp<<std::endl;
      int k=0;
      for(str=tmp;;str=NULL)
      {
        const char *token = strtok_r(str, ",", &saveptr);
        if(!token || !strlen(token)) break;
        if(k<limit)
          _s.push_back(dir+token);
        else
          _s.push_back(token);
        k++;
      }
      if(_s.empty()) continue;
      tbl.push_back(_s);
    }
  }
  
  void read_doubles_table(const char *fn,doubles_table& tbl)
  {
    tbl.clear();
    //read input file line by line, ignoring empty lines and lines starting with #
    std::ifstream in(fn);
    if(!in.good())
      REPORT_ERROR("Can't open file for reading!");
    
    while(in.good() && !in.eof())
    {
      char tmp[32000];
      in.getline(tmp,sizeof(tmp));
      if(!strlen(tmp)) continue;
      if(tmp[0]=='#') continue;
      char *str,*saveptr;
      doubles _s;
      //std::cout<<tmp<<std::endl;
      int k=0;
      for(str=tmp;;str=NULL)
      {
        const char *token = strtok_r(str, ",", &saveptr);
        if(!token || !strlen(token)) break;
        _s.push_back(atof(token));
        k++;
      }
      if(_s.empty()) continue;
      tbl.push_back(_s);
    }
  }
  
  bool table_verify(const doubles_table& tbl)
  {
    if(tbl.empty()) 
    {
      std::cerr<<"No data"<<std::endl;
      return false; //?
    }
    int sz=tbl[0].size();
    if(!sz)
    {
      std::cerr<<"First entry have no elements"<<std::endl;
      return false; //?
    }
    for(int i=0;i<tbl.size();++i)
    {
      if(tbl[i].size()!=sz)
      {
        std::cerr<<"Line :"<<i<<" have inconsistent number of entries:"
        <<tbl[i].size()<<" expected:"<<sz<<std::endl;
        for(int j=0;j<tbl[i].size();j++)
          std::cerr<<tbl[i][j]<<" - ";
        std::cerr<<std::endl;
        return false;
      }
    }
    return true;
  }
  
  //make sure all entries have same number of elements
  bool table_verify(const string_table& tbl)
  {
    if(tbl.empty()) 
    {
      std::cerr<<"No data"<<std::endl;
      return false; //?
    }
    int sz=tbl[0].size();
    if(!sz)
    {
      std::cerr<<"First entry have no elements"<<std::endl;
      return false; //?
    }
    for(int i=0;i<tbl.size();++i)
    {
      if(tbl[i].size()!=sz)
      {
        std::cerr<<"Line :"<<i<<" have inconsistent number of entries:"
            <<tbl[i].size()<<" expected:"<<sz<<std::endl;
        for(int j=0;j<tbl[i].size();j++)
          std::cerr<<tbl[i][j].c_str()<<" - ";
        std::cerr<<std::endl;
        return false;
      }
    }
    return true;
  }
  
  double variance(const minc::simple_volume<float>& v,const minc::simple_volume<float>& mean)
  {
    double r=0;
    for(int i=0;i<v.c_buf_size();i++)
    {
      double t=v.c_buf()[i]-mean.c_buf()[i];
      r+=t*t;
    }
    return r;
  }
  
  double sum(const minc_float_volume& v)
  {
    double r=0;
    for(int i=0;i<v.c_buf_size();i++)
    {
      double t=v.c_buf()[i];
      r+=t;
    }
    return r;
  }
  
  double sum(const minc_grid_volume& v)
  {
    double r=0;
    for(int i=0;i<v.c_buf_size();i++)
    {
      double t=v.c_buf()[i].sum();
      r+=t;
    }
    return r;
  }
  
  double sum(const minc_byte_volume& v)
  {
    double r=0;
    for(size_t i=0;i<v.c_buf_size();i++)
    {
      double t=v.c_buf()[i];
      r+=t;
    }
    return r;
  }
  
  double sum2(const minc_float_volume& v,const minc_byte_volume& mask)
  {
    if(v.size()!=mask.size())
      REPORT_ERROR("Volume size mismatch");
    
    double r=0;
    for(size_t i=0;i<v.c_buf_size();i++)
      if(mask.c_buf()[i])
    {
      double t=v.c_buf()[i];
      r+=t*t;
    }
    return r;
  }
  
  double sum2(const minc_float_volume& v)
  {
    double r=0;
    for(size_t i=0;i<v.c_buf_size();i++)
    {
      double t=v.c_buf()[i];
      r+=t*t;
    }
    return r;
  }
       
  double sum2(const minc_grid_volume& v,const minc_byte_volume& mask)
  {
    if(v.size()!=mask.size())
      REPORT_ERROR("Volume size mismatch");
    
    double r=0;
    for(size_t i=0;i<v.c_buf_size();i++)
      if(mask.c_buf()[i])
    {
      double t=v.c_buf()[i].mod2();
      r+=t;
    }
    return r;
  }

  void average_tables(string_table& tbl,volumes &v)
  {
    v.resize(tbl[0].size());
    minc_float_volume tmp;
    for(int i=0;i<tbl[0].size();i++) //go through modalities
    {
      minc_1_reader rdr1;
      rdr1.open(tbl[0][i].c_str());
      //rdr1.setup_read_float();
      load_simple_volume<float>(rdr1,v[i]);
      tmp.resize(v[i].size());
      for(int j=1;j<tbl.size();j++)
      {
        minc_1_reader rdr2;
        rdr2.open(tbl[j][i].c_str());
        load_simple_volume<float>(rdr2,tmp);
        v[i]+=tmp;
      }
      v[i]/=tbl.size();
    }
  }

  void average_tables(string_table& tbl,grids &v)
  {
    v.resize(tbl[0].size());
    minc_grid_volume tmp;
    for(int i=0;i<tbl[0].size();i++) //go through modalities
    {
      minc_1_reader rdr1;
      rdr1.open(tbl[0][i].c_str());
      //rdr1.setup_read_float();
      load_simple_volume(rdr1,v[i]);
      tmp.resize(v[i].size());
      for(int j=1;j<tbl.size();j++)
      {
        minc_1_reader rdr2;
        rdr2.open(tbl[j][i].c_str());
        load_simple_volume(rdr2,tmp);
        v[i]+=tmp;
      }
      v[i]/=IDX<float>(tbl.size(),tbl.size(),tbl.size());
    }
  }


  
  void mul(minc_float_volume &v2,const minc_float_volume& v1)
  {
    if(v1.size()!=v2.size())
      REPORT_ERROR("Volume size mismatch");
    
    for(int i=0;i<v2.c_buf_size();i++)
        v2.c_buf()[i]*=v1.c_buf()[i];
  }
  
  void masked_mul(minc_float_volume &v2,const minc_float_volume& v1,const minc_byte_volume &mask)
  {
    if(v1.size()!=v2.size())
      REPORT_ERROR("Volume size mismatch");
    if(v1.size()!=mask.size())
      REPORT_ERROR("Volume size mismatch");
    
    for(int i=0;i<v2.c_buf_size();i++)
      if(mask.c_buf()[i])
        v2.c_buf()[i]*=v1.c_buf()[i];
      else
        v2.c_buf()[i]=0;
  }


  
  void masked_mul(minc_grid_volume &v2,const minc_grid_volume& v1,const minc_byte_volume &mask)
  {
    if(v1.size()!=v2.size())
      REPORT_ERROR("Volume size mismatch");
    if(v1.size()!=mask.size())
      REPORT_ERROR("Volume size mismatch");
    
    for(int i=0;i<v2.c_buf_size();i++)
      if(mask.c_buf()[i])
        v2.c_buf()[i]*=v1.c_buf()[i];
      else
        v2.c_buf()[i]=IDX<float>(0,0,0);
  }

  
 
  std::string _dirname(const char *file)
  {
    char* tmp=new char[strlen(file)+1];
    strcpy(tmp,file);
    std::string r=::dirname(tmp);
    delete[] tmp;
    return r;
  }
  
  std::string _basename(const char *file)
  {
    char* tmp=new char[strlen(file)+1];
    strcpy(tmp,file);
    std::string r=::basename(tmp);
    delete[] tmp;
    return r;
  }

  int calc_selected(string_table& tbl,double threshold)
  {
    double total_sum=0.0;
  
    for(int i=1;i<tbl.size();i++)
    {
      total_sum+=fabs(atof(tbl[i][0].c_str()));
    }
    double sum=0.0;
    int selected;
    for(selected=1;selected<tbl.size();selected++)
    {
      sum+=fabs(atof(tbl[selected][0].c_str()));
      if((sum/total_sum)>threshold) break;
    }
    selected--;
    return selected;
  }

  void load_volume(const char *fn,minc::minc_grid_volume& vol)
  {
    minc_1_reader rdr;
    rdr.open(fn);
    load_simple_volume(rdr,vol);
  }
  
  void load_volume(const char *fn,minc::minc_byte_volume& vol)
  {
    minc_1_reader rdr;
    rdr.open(fn);
    load_simple_volume<unsigned char>(rdr,vol);
  }
  
  void load_volume(const char *fn,minc::minc_float_volume& vol)
  {
    minc_1_reader rdr;
    rdr.open(fn);
    load_simple_volume<float>(rdr,vol);
  }
  
  
  void load_volumes(strings& s,volumes& averages,int skip,bool verbose,bool ignore_missing)
  {
    averages.resize(s.size()-skip); //first column is eigenvalue
    for(int i=0;i<(s.size()-skip);i++)
    {
      minc_1_reader rdr;

      try
      {
        rdr.open(s[i+skip].c_str());
        load_simple_volume<float>(rdr,averages[i]);
      } catch(const minc::generic_error & err) {
        if(!ignore_missing) throw;
      }
    }
  }
  
  void check_volumes(volumes& vols,bool ignore_missing)
  {
    if(vols.empty()) return;
    
    if(ignore_missing) // go over all volumes and initialize missing ones with zeros
    {
      minc_float_volume::idx vol_size=IDX<size_t>(0,0,0);
      for(int i=0;i<vols.size();i++)
      {
        if(!vols[i].empty())
        {
          vol_size=vols[i].size();
          break;
        }
      }
      
      if(!vol_size[0]) //all volumes are empty?
        return;
      for(int i=0;i<vols.size();i++)
      {
        if(vols[i].empty())
        {
          vols[i].resize(vol_size);
          vols[i]=0.0;//fill with zeroes
        }
      }
    }
    
    for(int i=1;i<vols.size();i++)
    {
      if(vols[i].size()!=vols[0].size())
        REPORT_ERROR("Volumes size mismatch!");
    }
  }

  void load_all_volumes(string_table & tbl,std::vector<volumes>& all_volumes,int selected,bool verbose)
  {
    all_volumes.resize(selected);
    for(int i=0;i<selected;i++)
    {
      all_volumes[i].resize(tbl[0].size()-1);
      load_volumes(tbl[i+1],all_volumes[i],1,verbose);
    }
  }
  
  void load_all_volumes_0(string_table & tbl,std::vector<volumes>& all_volumes,int selected,bool verbose)
  {
    all_volumes.resize(selected);
    for(int i=0;i<selected;i++)
    {
      all_volumes[i].resize(tbl[0].size());
      load_volumes(tbl[i],all_volumes[i],0,verbose);
    }
  }

  void substract(volumes& v1,const volumes& v2)
  {
    if(v2.size()<v1.size())
      REPORT_ERROR("Too many input volumes");
  
    for(int i=0;i<v1.size();i++)
      v1[i]-=v2[i];
  }


  void save_volume(minc::minc_grid_volume& vol,const char *fn,const char *like,const std::string & append_history)
  {
    minc_1_reader rdr;
    rdr.open(like,false,true);// I need only the metadata
    minc_1_writer wrt;
    wrt.open(fn,rdr.info(),2,NC_FLOAT);
    
    wrt.copy_headers(rdr);
    if(! append_history.empty())
      wrt.append_history(append_history.c_str());
    
    save_simple_volume(wrt,vol);
  }
  
  void save_volume(minc::minc_float_volume& vol,const char *fn,const char *like,const std::string & append_history)
  {
    minc_1_reader rdr;
    rdr.open(like,false,true);// I need only the metadata
    minc_1_writer wrt;
    wrt.open(fn,rdr.info(),2,NC_FLOAT);
    
    wrt.copy_headers(rdr);
    if(! append_history.empty())
      wrt.append_history(append_history.c_str());
    
    save_simple_volume<float>(wrt,vol);
  }
  
  void save_volume(minc::minc_byte_volume& vol,const char *fn,const char *like,const std::string & append_history)
  {
    minc_1_reader rdr;
    rdr.open(like,false,true);// I need only the metadata
    minc_1_writer wrt;
    wrt.open(fn,rdr.info(),2,NC_BYTE,false);
    
    wrt.copy_headers(rdr);
    if(! append_history.empty())
      wrt.append_history(append_history.c_str());
    
    save_simple_volume<unsigned char>(wrt,vol);
  }
  
  void save_volumes(strings& s,volumes& vols,const char *like,const std::string & append_history)
  {
    for(int i=0;i<s.size();i++)
    {
      save_volume(vols[i],s[i].c_str(),like,append_history);
    }
  }

  void load_volumes(strings& s,grids& vols,int skip,bool verbose)
  {
    vols.resize(s.size()-skip); //first column is eigenvalue
    for(int i=0;i<(s.size()-skip);i++)
    {
      /*if(verbose) 
        std::cout<<s[i+skip].c_str()<<" "<<std::flush;*/
      
      minc_1_reader rdr;
      rdr.open(s[i+skip].c_str());
      load_simple_volume(rdr,vols[i]);
    }
  }

  void load_all_volumes(string_table & tbl,std::vector<grids>& all_volumes,int selected)
  {
    all_volumes.resize(selected);
    for(int i=0;i<selected;i++)
    {
      all_volumes[i].resize(tbl[0].size()-1);
      load_volumes(tbl[i+1],all_volumes[i]);
    }
  }

  void substract(grids& v1,const grids& v2)
  {
    if(v2.size()<v1.size())
      REPORT_ERROR("Too many input volumes");
  
    for(int i=0;i<v1.size();i++)
      v1[i]-=v2[i];
  }
  
  double calc_dv(const char* minc_file)
  {
    minc_1_reader rdr;
    rdr.open(minc_file,false,true);// I need only the metadata
    return rdr.nspacing(1)*rdr.nspacing(2)*rdr.nspacing(3);
  }

  //based on the example from libstdc
  inline bool
      isws (char c, char const * const wstr)
  {
    return (strchr(wstr,c) != NULL);
  }
  
/*****************************************************************
  * Simplistic and quite Standard, but a bit slow.  This should be
  * templatized on basic_string instead, or on a more generic StringT
  * that just happens to support ::size_type, .substr(), and so on.
  * I had hoped that "whitespace" would be a trait, but it isn't, so
  * the user must supply it.  Enh, this lets them break up strings on
  * different things easier than traits would anyhow.
   void stringtok (Container &l, std::string const &s, char const * const ws = " \t\n")

*/
  void string_to_vector(const std::string& s,std::vector<double>& l,const char *ws/*=" ,"*/)
  {
    const std::string::size_type  S = s.size();
    std::string::size_type  i = 0;

    while (i < S) {
      // eat leading whitespace
      while ((i < S) && (isws(s[i],ws)))  ++i;
      if (i == S)  return;  // nothing left but WS

      // find end of word
      std::string::size_type  j = i+1;
      while ((j < S) && (!isws(s[j],ws)))  ++j;

      // add word
      l.push_back(atof(s.substr(i,j-i).c_str()));

        // set up for next loop
      i = j+1;
    }
  }
  
  double volume_rms_diff(const minc::minc_float_volume& in1,const minc::minc_float_volume& in2)
  {
    if(in1.size()!=in2.size())
      REPORT_ERROR("Size mismatch!");
    double diff=0;
    int cnt=0;
  
    for(int i=0;i<in1.c_buf_size();i++)
    {
      double d=in1.c_buf()[i]-in2.c_buf()[i];
      diff+=d*d;
    //w+=g2.c_buf()[i]*g2.c_buf()[i];
      cnt++;
    }
    if(cnt>0) diff/=cnt;
    return sqrt(diff);
  }  
  
}; //minc
