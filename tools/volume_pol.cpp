#include "minc_wrappers.h"
#include "minc_helpers.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>
//#include "data_proc.h"

#include "histograms.h"
#include "gsl_glue.h"
#include "gsl_gauss.h"

namespace minc
{
	template < class T > void build_histogram(histogram< T > &hist,minc::image3d::Pointer img)
	{
	
		//1. get min, max
		voxel_type min, max;
		int count = get_image_limits(img,min,max);
	
		hist.clear();
		hist.set_limits(min, max);
		minc::image3d_iterator it(img, img->GetLargestPossibleRegion ());
    int cnt=0;
		//2. populate histogram
		for (it.GoToBegin (); !it.IsAtEnd (); ++it) {
			T v=it.Value();
			if (fabs(v) > 1e20) //remove minc undefined voxels
				continue;
			hist[v]++;
      cnt++;
		}
		hist /= cnt;
	}
	
	template < class T > void build_histogram (histogram < T > &hist,
												minc::image3d::Pointer img,
												minc::mask3d::Pointer msk)
	{
		//1. get min, max
		minc::image3d_iterator it (img, img->GetLargestPossibleRegion ());
		minc::mask3d_const_iterator it3 (msk, msk->GetLargestPossibleRegion ());
		int count = 0;
	
		T 	min = std::numeric_limits < T >::max (),
				max = -std::numeric_limits < T >::max ();
		it3.GoToBegin();
		for (it.GoToBegin (); !it.IsAtEnd (); ++it) {
			T v=it.Value();
			if (fabs(v) > 1e20)
				continue;
			if (it3.Value()) {
				count++;
				if (v < min)
					min = v;
				if (v > max)
					max = v;
			}
			++it3;
		}
		hist.clear ();
		hist.set_limits (min, max);
		//2. populate histogram
		it3.GoToBegin ();
		for (it.GoToBegin (); !it.IsAtEnd (); ++it) {
			T v=it.Value ();
			if (fabs(v) > 1e20)
				continue;
	
			if (it3.Value ())
				hist[v]++;
			++it3;
		}
		hist /= count;
	}
	
	template < class T > void build_joint_histogram (joint_histogram < T > &hist,
															minc::image3d::Pointer img1,
															minc::image3d::Pointer img2)
	{
		if(img1->GetLargestPossibleRegion ().GetSize()!=
						img2->GetLargestPossibleRegion ().GetSize())
			REPORT_ERROR("Volume dimensions mismatch");
	
		//1. get min, max
		minc::image3d_iterator it1 (img1, img1->GetLargestPossibleRegion ());
		minc::image3d_iterator it2 (img2, img2->GetLargestPossibleRegion ());
	
		//    value_type min1,min2,max1,max2;
	
		//    int count = get_image_limits(img1,min1,max1);
		//    get_image_limits(img2,min2,max2);
		int count=0;
	
		//hist.set_limits(min1,max1);
		//hist.set_joint_limits(min2,max2);
		hist=T(0);
	
		//2. populate histogram
		for (it1.GoToBegin (),it2.GoToBegin(); !it1.IsAtEnd (); ++it1,++it2) {
			T v1,v2;
			if ((v1=it1.Value ()) > 1e10 || (v2=it2.Value())>1.e10)
				continue;
			hist[v1]++;
			hist(v1,v2)++;
			count++;
		}
		hist.normalize(count);
	}
};

void print_expression (std::ostream & out, const std::vector < double >&coeff2,
                  float min2, float max2,bool noclamp=false)
{
	int order = coeff2.size ();
	out << "clamp(";
  switch(order)
  {
    case 1:
      out << coeff2[0]<<"*A[0]";break;
    case 2:
      out << coeff2[0]<<"*A[0]+"<<coeff2[1];break;
    default:
    for (int i = 0; i < order; i++)
    {
      out << coeff2[i];
      for (int j = 0; j < i; j++)
        out << "*A[0]";
      if (i != (order - 1))
        out << "+";
    }
  }
	out << "," << min2 << "," << (noclamp?max2*1000:max2) << ")" << std::endl;
	out.flush();
}


using namespace  std;
using namespace  minc;
const int buckets = 100;

void show_usage(void)
{
	std::cerr << "Usage: volume_pol <input> <template> " << endl
	<< "[--verbose] be verbose" << endl
	<< "[--quiet] be quiet " << endl
	<< "[--clobber] clobber the output files" << endl
	<< "[--source_mask <mask>] use mask for source file" << endl
  << "[--target_mask <mask>] use mask for tamplate file" << endl
	<< "[--order <n>] approximation order" << endl 
  << "[--hist <histogram.dat>] output histograms" << endl
	<< "[--joint <joint histogram>] produce joint histogram (files should have same dimensions)"<<endl
	<< "[--expfile <minccalc expression file>] write output to the expression file " << endl
	<< "[--version] print version" << endl
  << "[--min <min>] map 0.1% to min" << endl
  << "[--max <max>] map 99.9% to max" << endl
  << "[--noclamp] don't clamp highlights" <<endl;
}


typedef double histogram_type;

int main (int argc, char **argv)
{
	histogram<histogram_type>
    hist1(buckets),
    hist2(buckets),
    hist3(buckets),
    hist4(buckets);
    
	std::vector<histogram<histogram_type> >  joint(buckets, histogram < histogram_type >(buckets));
	int verbose = 0;
	int clobber = 0;
	int histogram = 1;
  int noclamp=0;
  histogram_type map_min=0,map_max=0;
	bool use_mask = false, make_output_file = false;
	int  order = 4;
	std::string source_mask_file,target_mask_file;
	std::string hist_file, joint_hist_file;
	std::string outfile, infile, templatefile, exp_file;
  
	static struct option long_options[] =
	  {
		  {"verbose", no_argument, &verbose, 1},
		  {"quiet", no_argument, &verbose, 0},
		  {"clobber", no_argument, &clobber, 1},
		  {"source_mask", required_argument, 0, 'm'},
      {"target_mask", required_argument, 0, 't'},
		  {"order", required_argument, 0, 'o'},
		  {"hist", required_argument, 0, 'h'},
		  {"joint",     required_argument,0, 'j'},
		  {"expfile", required_argument, 0, 'e'},
		  {"version", no_argument, 0, 'v'},
      {"min", required_argument, 0, 'i'},
      {"max", required_argument, 0, 'a'},
      {"noclamp",no_argument,&noclamp,1},
		  {0, 0, 0, 0}
	  };

	int c;
	for (;;)
	{
		/* getopt_long stores the option index here. */
		int
		option_index = 0;

		c = getopt_long (argc, argv, "vqm:o:vh:j:e:i:a:t:n", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
		case 'm':
			source_mask_file = optarg;
			break;
		case 't':
			target_mask_file = optarg;
			break;
		case 'o':
			order = atoi(optarg);
			break;
		case 'h':
			hist_file = optarg;
			break;
		case 'j':
			joint_hist_file=optarg;
			break;
		case 'e':
			exp_file = optarg;
			break;
		case 'v':
			cout << "Version: 1.1" << endl;
			return 0;
    case 'i':
      map_min=atof(optarg);
      break;
    case 'a':
      map_max=atof(optarg);
      break;
		case '?':
			/* getopt_long already printed an error message. */
		default:
			show_usage ();
			return 1;
		}
	}

	if ((argc - optind) < 2)
	{
		show_usage ();
		return 1;
	}
	infile = argv[optind];
	templatefile = argv[optind + 1];
	if ((argc - optind) > 2)
	{
		outfile = argv[optind + 2];
		make_output_file = true;
	}

	try
	{
		image3d::Pointer img1   (image3d::New()), img2(image3d::New());
		mask3d::Pointer  src_msk(mask3d::New()),  trg_msk(mask3d::New());
		if (order <= 0)
			order = 1;

    load_minc(infile.c_str (), img1);
    
		image3d::SizeType s1 = img1->GetLargestPossibleRegion().GetSize ();
		if(verbose)
			cout<<"Image "<<s1[0]<<" x "<<s1[1]<< " x "<<s1[2]<<endl;

		load_minc(templatefile.c_str (), img2);

    image3d::SizeType s2 = img2->GetLargestPossibleRegion().GetSize ();
		if(verbose)
			cout<<"Template "<<s2[0]<<" x "<<s2[1]<< " x "<<s2[2]<<endl;

		if(!joint_hist_file.empty() && s1!=s2)
		{
			cerr<<"For Joint Histogram Images should have the same dimensions, use mincresample!"<<endl;
			return 1;
		}

		if (!source_mask_file.empty())
		{
  		load_minc(source_mask_file.c_str (), src_msk);
			mask3d::SizeType ms = src_msk->GetLargestPossibleRegion().GetSize ();
			if (ms != s1)
			{
				cerr << "Mask is wrong size, use mincresample!" << endl;
				return 1;
			}
		}
    
		if (!target_mask_file.empty())
		{
  		load_minc(target_mask_file.c_str (), trg_msk);
			mask3d::SizeType ms  = trg_msk->GetLargestPossibleRegion().GetSize ();
			if (ms != s2)
			{
				cerr << "Mask is wrong size, use mincresample!" << endl;
				return 1;
			}
		}

		MNK_Gauss_Polinomial_Mod pol(order);

		if (!source_mask_file.empty())
			build_histogram<histogram_type>(hist1, img1, src_msk);
		else
			build_histogram<histogram_type>(hist1, img1);

		if (verbose)
			cout << "Image min=" << hist1.min () << " max=" << hist1.max() << endl;

		if (!target_mask_file.empty())
			build_histogram<histogram_type>(hist2, img2, trg_msk);
		else
			build_histogram<histogram_type>(hist2, img2);

		if (verbose)
			cout << "Template min=" << hist2.min() << " max=" << hist2.max() << endl;

		//MNK_Gauss_Polinomial pol2(order);
    MNK_Gauss_Polinomial pol2(order); 


    histogram_type omin,omax;
    //map ramge
    if(map_min<map_max)
    {
      omin=map_min;
      omax=map_max;
      histogram_type tmin=hist2.find_percentile(0.001);
      histogram_type tmax=hist2.find_percentile(0.999);
      
      histogram_type k=(map_max-map_min)/(tmax-tmin);
      
      for (histogram_type pc = 0.001; pc < 1.0; pc += 0.001)
        pol.accumulate(hist1.find_percentile(pc), k*(hist2.find_percentile(pc)-tmin)+map_min);

    } else {
      omin=hist2.min ();
      omax=hist2.max ();
      for (histogram_type pc = 0.001; pc < 1.0; pc += 0.001)
      {
        pol.accumulate(hist1.find_percentile(pc), hist2.find_percentile(pc));
      }
    }
		std::vector < double > coeff2 (order);

    pol.solve(coeff2);
    //if(verbose)
		// cout << "Condition number:"<<cond << endl;

		if (verbose)
		{
			cout << "Histogram fitting:" << endl;
			for (int i = 0; i < order; i++)
				cout << coeff2[i] << ",";

			//cout << "rank=" << rank2 << endl;
		}
		if (exp_file.empty ())
		{
			if (verbose)
				cout << "minccalc string:" << endl;
			print_expression (cout, coeff2, omin, omax,noclamp);
		}
		else
		{
			if (!clobber && !access (exp_file.c_str (), F_OK))
			{
				cerr << exp_file.c_str () << " Exists!" << endl;
				return 1;
			}
			if (verbose)
				cout << "Writing output to " << exp_file.c_str () << endl;
			ofstream
			o_exp (exp_file.c_str ());
			print_expression (o_exp, coeff2, omin, omax,noclamp);
		}

		if (!hist_file.empty ())
		{
			if (!clobber && !access (hist_file.c_str (), F_OK))
			{
				cerr << hist_file.c_str () << " Exists!" << endl;
				return 1;
			}
			ofstream
			out1 (hist_file.c_str ());
			if (out1.bad ())
			{
				cerr << "Error writing to:" << hist_file.c_str () << endl;
				return 1;
			}
      
			out1 << "#value,hist_orig,prob_lsq_fit" << endl;
			for (int i = 0; i < buckets; i++)
			{
				out1<< hist1.value(i) << " "<< hist1[i]<< " "<< /*hist4[i]<< " "<<*/ pol.fit(coeff2, hist1.value(i), hist2.min(), hist2.max())<< endl;
				//out_joint<<hist1.value(i)<<" ";
			}
		}


		if(!joint_hist_file.empty())
		{
			if(!clobber && !access(joint_hist_file.c_str(),F_OK))
			{
				cerr<<joint_hist_file.c_str()<<" Exists!"<<endl;
				return 1;
			}

			joint_histogram<histogram_type> j_hist(buckets,hist1.min(),hist1.max());
			j_hist.set_joint_limits(hist2.min(),hist2.max());
			build_joint_histogram(j_hist,img1,img2);

			ofstream out_joint(joint_hist_file.c_str());
			if(out_joint.bad())
			{
				cerr<<"Error writing to:"<<joint_hist_file.c_str()<<endl;
				return 1;
			}
			j_hist.save(out_joint);

		}
	}
  
	catch (const minc::generic_error & err)
	{
		cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
		return 1;
	}
}
