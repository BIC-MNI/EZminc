#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <math.h>
#include <limits>
#include <unistd.h>
//#include "data_proc.h"

#include <minc_histograms.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>

#include "gsl_glue.h"
#include "gsl_gauss.h"
#include <time_stamp.h>    // for creating minc style history entry

void print_expression (std::ostream & out, const std::vector < double >&coeff2,
                  float min2, float max2,bool noclamp=false)
{
	int order = coeff2.size ();
  if(!noclamp)
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
  if(!noclamp)
    out << "," << min2 << "," << (noclamp?max2*1000:max2) << ")";
  out << std::endl;
  
	out.flush();
}


using namespace  std;
using namespace  minc;
const int buckets = 100;

void show_usage(void)
{
	std::cerr << "Usage: volume_pol <input> <template> [output_file]" << endl
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
  << "[--noclamp] don't clamp highlights" <<endl
  << "[--kl] calculate Kullback–Leibler divergence instead of fitting"<<endl
  << "[--ks] calculate Kolmogorov–Smirnov distance instead of fitting"<<endl;
}


typedef double histogram_type;

int main (int argc, char **argv)
{
    
  
	int verbose = 0;
	int clobber = 0;
	int histogram = 1;
  int noclamp=0;
  histogram_type map_min=0,map_max=0;
	bool use_mask = false;
	int  order = 4;
  int calc_kl=0,calc_ks=0;
	std::string source_mask_file,target_mask_file;
	std::string hist_file, joint_hist_file;
	std::string outfile, infile, templatefile, exp_file;
  
  char *_history = time_stamp(argc, argv); 
  std::string history=_history;
  
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
      {"kl",no_argument,&calc_kl,1},
      {"ks",no_argument,&calc_ks,1},
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
    
    if (!clobber && !access (outfile.c_str (), F_OK))
    {
      cerr << outfile.c_str () << " Exists!" << endl;
      return 1;
    }
  }
  
  if (!exp_file.empty() && !clobber && !access (exp_file.c_str (), F_OK))
  {
    cerr << exp_file.c_str () << " Exists!" << endl;
    return 1;
  }
  
  if(!joint_hist_file.empty() && !clobber && !access(joint_hist_file.c_str(),F_OK))
  {
    cerr<<joint_hist_file.c_str()<<" Exists!"<<endl;
    return 1;
  }
  

	try
	{
    minc::histogram<histogram_type>  hist1(buckets),hist2(buckets);
        
    simple_volume<double> img1,img2;
    minc_byte_volume src_msk,trg_msk;

		
    if (order <= 0)
			order = 1;

    minc_1_reader rdr1;
    rdr1.open(infile.c_str());
    load_simple_volume<double>(rdr1,img1);
    

    minc_1_reader rdr2;
    rdr2.open(templatefile.c_str());
    load_simple_volume<double>(rdr2,img2);

		if(!joint_hist_file.empty() && img1.size()!=img2.size())
		{
			cerr<<"For Joint Histogram Images should have the same dimensions, use mincresample!"<<endl;
			return 1;
		}

		if (!source_mask_file.empty())
		{
      minc_byte_volume src_mask;
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<source_mask_file.c_str()<<std::endl;
      
      rdr.open(source_mask_file.c_str());
      load_simple_volume(rdr,src_msk);
      
      if(src_msk.size()!=img1.size())
        REPORT_ERROR("Source mask size mismatch");
    }
    
		if (!target_mask_file.empty())
		{
      minc_byte_volume src_mask;
      minc_1_reader rdr;
      if(verbose) std::cout<<"loading mask:"<<target_mask_file.c_str()<<std::endl;
      
      rdr.open(target_mask_file.c_str());
      load_simple_volume(rdr,trg_msk);
      
      if(trg_msk.size()!=img2.size())
        REPORT_ERROR("Target mask size mismatch");
		}
		
		int size_hist1=0,size_hist2=0;

		if (!source_mask_file.empty())
      size_hist1=build_histogram(hist1,img1,src_msk);
		else
      size_hist1=build_histogram(hist1,img1);

		if (verbose)
			cout << "Image min=" << hist1.min () << " max=" << hist1.max() << endl;

		if (!target_mask_file.empty())
      size_hist2=build_histogram(hist2,img2,trg_msk);
		else
      size_hist2=build_histogram(hist2,img2);

		if (verbose)
			cout << "Template min=" << hist2.min() << " max=" << hist2.max() << endl;

    MNK_Gauss_Polinomial_Mod pol(order);
    MNK_Gauss_Polinomial pol2(order); 
    
    if(calc_kl)
    {
      if(verbose)
        std::cout<<"Kullback–Leibler divergence:";
      std::cout<<kl_distance(hist1,hist2)<<std::endl;
      return 0;
    } else if(calc_ks) {
      if(verbose)
        std::cout<<"Kolmogorov–Smirnov distance:";
      double dist=ks_distance(hist1,hist2);
      std::cout<<dist<<std::endl;
/*      if(verbose)
        std::cout<<"p="<<ks_significance(dist,size_hist1,size_hist2)<<std::endl;*/
      return 0;
    }


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
			minc::joint_histogram<histogram_type> j_hist(buckets,hist1.min(),hist1.max());
			j_hist.set_joint_limits(hist2.min(),hist2.max());
      
      build_joint_histogram<double>(j_hist,img1,img2);

			ofstream out_joint(joint_hist_file.c_str());
			if(out_joint.bad())
			{
				cerr<<"Error writing to:"<<joint_hist_file.c_str()<<endl;
				return 1;
			}
			j_hist.save(out_joint);
		}
    
    if(!outfile.empty())
    {
      //convert values 
      for(int i=0;i<img1.c_buf_size();i++)
      {
        if(noclamp)
          img1.c_buf()[i]=pol.fit(coeff2, img1.c_buf()[i]);
        else
          img1.c_buf()[i]=pol.fit(coeff2, img1.c_buf()[i], hist2.min(), hist2.max());
      }
      
      if (verbose)
        cout << "Writing output to " << outfile.c_str () << endl;
      
      minc_1_writer wrt;
    
      wrt.open(outfile.c_str(),rdr1);
      wrt.append_history(history.c_str());
    
      save_simple_volume<double>(wrt,img1);
      
    }
	}
  
	catch (const minc::generic_error & err)
	{
		cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
		return 1;
	}
}
