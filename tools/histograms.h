#ifndef __HISTOGRAMS_H__
#define __HISTOGRAMS_H__

#include <valarray>
#include <vector>

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
			return T (i) / _k + _min;
		}
	
		histogram (int buckets, T min = 0, T max = 1):
				_hist (buckets), _min (min), _max (max),
				_range (_max - _min)
		{
			_hist = 0.0;
			_k = T (_hist.size ()) / _range;
			;
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
			_min = a._min;
			_max = a._max;
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
			int i, j;
			for (i = 0, j = 0; i < size(); i++) {
				acc += (*this)[i];
				if (acc >= pc)
					break;
				if ((*this)[i] > 0.0)
					j = i;
			}
			if(i > 0)
      {
        T k=(acc-pc)/(*this)[i];
        return value(i)*(1.0-k) + value(j)*k;
        //return value(j);
      }
			else
				return value(0);
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
	
	
}; //minc


#endif //__HISTOGRAMS_H__
