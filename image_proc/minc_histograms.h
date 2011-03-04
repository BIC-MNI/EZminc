/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: image histogram manipulation routines
@COPYRIGHT  :
              Copyright 2009 Vladimir Fonov, 
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
#ifndef __MINC_HISTOGRAMS_H__
#define __MINC_HISTOGRAMS_H__

#include <valarray>
#include <vector>
#include <minc_io_simple_volume.h>
#include <iostream>
#include <fstream>

namespace minc
{
  
	template<class T> class histogram
	{
	protected:
		std::valarray< T > _hist;
		T _min, _max, _range;
		T _k;
	
		int _idx (T idx) const
		{
			int i = int (floor (_k * (idx - _min)));
      
			if (i < 0)
				i = 0;
			if (i >= _hist.size ())
				i = _hist.size () - 1;
			return i;
		}
	
		int _idx (int i) const
		{
			if (i < 0)
				i = 0;
			if (i >= _hist.size ())
				i = _hist.size () - 1;
			return i;
		}
	public:
	
		int idx (T val) const
		{
			return _idx (val);
		}
	
		T value (int i) const
		{
			return T(i)/_k + _min + 0.5/_k;
		}
	
		histogram (int buckets, T min = 0, T max = 1):
				_hist (buckets), _min (min), _max (max),
				_range (_max - _min)
		{
			_hist = 0.0;
			_k = T (_hist.size ()) / _range;
  	}
	
		void set_limits(T min, T max)
		{
			_max = max;
			_min = min;
			_range = _max - _min;
			_k = T (_hist.size ()) / _range;
		}
	
		T range() const
		{
			return _range;
		}
	
		T min(void) const
		{
			return _min;
		}
	
		T max(void) const
		{
			return _max;
		}
	
		T& operator[](T idx)
		{
			return _hist[_idx(idx)];
		}
	
		T operator[](T idx) const
		{
			return _hist[_idx (idx)];
		}
	
		T & operator[](int i)
		{
			return _hist[_idx (i)];
		}
	
		T operator[](int i) const
		{
			return _hist[_idx (i)];
		}
	
		const histogram & operator=(histogram & a)
		{
      set_limits(a._min,a._max);
			_hist = a._hist;
			return *this;
		}
	
		const histogram & operator= (T v)
		{
			_hist = v;
			return *this;
		}
	
		const histogram & operator/= (T v)
		{
			_hist /= v;
			return *this;
		}
	
		const histogram & operator*= (T v)
		{
			_hist *= v;
			return *this;
		}
	
		int size (void) const
		{
			return _hist.size ();
		}
	
		T find_percentile(T pc) const
		{
			T acc = 0.0;
			int i, j, k;
      
			for (i = 0, j = 0; i < size(); i++) 
      {
        T prob=(*this)[i];
				acc += prob;
        
        //std::cout<<"\tCount:"<<i<<" prob="<<acc*100<<std::endl;
        
				if (acc >= pc)
					break;
          
				if (prob > 0.0)
					j = i;
			}
      
			if(i > 0 ) //normal case
      {
        T k=(acc-pc)/(*this)[i];
        
        return value(i)*(1.0-k) + value(j)*k+0.5/_k;
        //return value(i);
        
      } else {
        T k=(acc-pc)/(*this)[0];
        
        return (value(0)+0.5/_k)*(1.0-k) + _min*k;
        //return value(0);
      }
		}
    
    T bimodalT(void)
    {
    //using Otsu'97 algorithm
      int thr=0;
      T mua,sa;
      T mub,sb;
      T min_s=1e10;
      T cnta,cntb;
      T s;
      
      for (int i = 1; i < (size()-1); i++) {
        mua=mub=0;
        sa=sb=0;
        cnta=cntb=0;
        s=0;
        
        for (int j = 0; j < i; j++)
        {
          mua+=value(j)*_hist[j];
          cnta+=_hist[j];
        }
        if(cnta>0)
          mua/=cnta;
          
        for (int j = 0; j < i; j++)
          sa+=(value(j)-mua)*(value(j)-mua)*_hist[j];
          
        if(cnta>0)
          sa/=cnta;
        
        for (int j = i; j < size(); j++)
        {
          mub+=value(j)*_hist[j];
          cntb+=_hist[j];
        }
        
        if(cntb>0)
          mub/=cntb;
          
        for (int j = i; j < size(); j++)
          sb+=(value(j)-mub)*(value(j)-mub)*_hist[j];
        sb/=cntb;
        
        s=(cntb*sb+cnta*sa)/(cnta+cntb);
        
        if(s<min_s)
        {
          min_s=s;
          thr=i;
        }
        }
      return value(thr);
    }
	
		void clear (void)
		{
			_hist = 0;
		}
	
		void save(std::ostream& out) const
		{
			for(int i=0;i<size();i++)
				out<<value(i)<<" "<<_hist[i]<<std::endl;
		}
    
    void save(const char *file) const
    {
      std::ofstream out(file);
      save(out);
    }
            
    void seed(const T& val,int count=1)
    {
      int i = idx(val);
      
      if (i >= (_hist.size()-1)) {
        i = _hist.size () - 1;
        _hist[i]+=count;
      } else {
        
        T d1=fabs(value(i)-val);
        T d2=fabs(value(i+1)-val);
        
        _hist[i]+=d2*count/(d1+d2);
        _hist[i+1]+=d1*count/(d1+d2);
        //_hist[i]+=count;
      }
    }

	};

	template< class T> class joint_histogram: public histogram<T>
	{
		typedef histogram<T> _histogram;
		typedef std::vector< histogram<T> > _histograms;
	protected:
		_histograms _joint;
		class _apply_min_max
		{
			const T &_min;
			const T &_max;
		public:
	
			_apply_min_max(const T& min,const T& max):_min(min),_max(max)
			{}
	
			void operator()(_histogram& h)
			{
				h.set_limits(_min,_max);
			}
		};
	
		class _apply_norm
		{
			const T &_norm;
		public:
			_apply_norm(const T& norm):_norm(norm)
			{ }
	
			void operator()(_histogram& h)
			{
				h/=_norm;
			}
		};
		class _out_file
		{
			std::ostream &_out;
		public:
	
			_out_file(std::ostream &out):_out(out)
			{}
	
			void operator()(const _histogram& h)
			{
				for(int i=0;i<h.size();i++)
					_out<<h[i]<<" ";
				_out<<std::endl;
			}
		};
		class _set_value
		{
			const T& _v;
		public:
			_set_value(const T& v):_v(v)
			{}
			void operator()(_histogram& h)
			{
				h=_v;
			}
		};
	public:
		joint_histogram(int buckets, T min = 0, T max = 1):
				histogram<T>(buckets,min,max),
				_joint(buckets,histogram<T>(buckets,min,max))
		{}
	
		void set_joint_limits(T min,T max)
		{
			_apply_min_max op(min,max);
			std::for_each(_joint.begin(),_joint.end(),op);
		}
	
		T& operator()(T x, T y)
		{
			return _joint[idx(x)][y];
		}
	
		const T& operator()(T x, T y) const
		{
			return _joint[idx(x)][y];
		}
	
		T& operator()(int x, int y)
		{
			return _joint[x][y];
		}
	
		const T& operator()(int x, int y) const
		{
			return _joint[x][y];
		}
	
		void normalize(T norm)
		{
			_apply_norm op(norm);
			std::for_each(_joint.begin(),_joint.end(),op);
		}
	
		_histogram& hist(T v)
		{
			return _joint[idx(v)];
		}
	
		const _histogram& hist(T v) const
		{
			return _joint[idx(v)];
		}
	
		_histogram& hist(int v)
		{
			return _joint[v];
		}
	
		const _histogram& hist(int v) const
		{
			return _joint[v];
		}
	
		void save(std::ostream &out) const
		{
			_out_file op(out);
			std::for_each(_joint.begin(),_joint.end(),op);
		}
	
		const joint_histogram<T>& operator=(T v)
		{
			_set_value op(v);
			std::for_each(_joint.begin(),_joint.end(),op);
		}
	
	};
	
  template < class T,class D > void build_histogram(histogram<D > &hist,const simple_volume<T> & img)
  {
		//1. get min, max
    T img_min, img_max;
    volume_min_max(img,img_min,img_max);
    
    int count = img.c_buf_size();
	
    hist.clear();
    hist.set_limits(img_min, img_max);
    int cnt=0;
		//2. populate histogram
    for (int i=0; i<count; i++) {
      hist.seed(img.c_buf()[i]);
      cnt++;
    }
    hist /= cnt;
  }
	
  template < class T,class D> void build_histogram (histogram < D > &hist,
                                             const simple_volume<T> & img,
                                             const simple_volume<unsigned char> & mask)
  {
	
    T min,max;
    volume_min_max(img,mask,min,max);
    
    hist.clear();
    hist.set_limits (min, max);
		//2. populate histogram
    
    int cnt=0;
    for (int i=0; i<img.c_buf_size(); i++) {
 
      if (mask.c_buf()[i])
      {
        hist.seed(img.c_buf()[i]);
        cnt++;
      }
    }
    hist /= cnt;
  }
	
  template < class T > void build_joint_histogram (joint_histogram < T > &hist,
      const simple_volume<T> img1,
      const simple_volume<T> img2)
  {
    if(img1.size()!=img2.size())
      REPORT_ERROR("Volume dimensions mismatch");
	
    T img1_min,img1_max;
    T img2_min,img2_max;
    
    volume_min_max(img1,img1_min,img1_max);
    volume_min_max(img2,img2_min,img2_max);
    
    int count=0;
	
		hist.set_limits(img1_min,img1_max);
		hist.set_joint_limits(img2_min,img2_max);
    
    hist=T(0);
	
		//2. populate histogram
    for (int i=0; i<img1.c_buf_size(); i++) {
      hist[img1.c_buf()[i]]++;
      hist(img1.c_buf()[i],img2.c_buf()[i])++;
      count++;
    }
    hist.normalize(count);
  }
  
  //! estimate sample mu using discrete classes
  void estimate_mu(const histogram<double>& input,
                 const std::vector<int>& cls,
                 std::vector<double>&  mu );
                 
  //!  apply k-means classify to histogram
  void apply_k_means_classify(const histogram<double>& input,
                      const std::vector<double>&  mu,
                      std::vector<int>& cls );
                      
  //!  apply k-means classify to volume
  void apply_hard_classify(simple_volume<float> & input,
                            minc_byte_volume    &mask,
                            std::vector<double>  mu,
                            minc_byte_volume& cls
                          );
                           
  //! perform simple k-means classify
  void simple_k_means( histogram<double>& hist, std::vector<double> & mu,int k_means,int maxiter);
  

}; //minc
#endif //__MINC_HISTOGRAMS_H__
