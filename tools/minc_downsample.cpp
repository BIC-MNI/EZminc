#include <iostream>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>

#include <time_stamp.h>    // for creating minc style history entry
// for get_opt_long
#include <getopt.h>
#include <cstdlib>
#include <unistd.h>

using namespace  std;
using namespace  minc;


void show_usage(const char *name)
{
  std::cerr 
	  << "Usage: "<<name<<" <input> <output> " << std::endl
    << "Optional parameters:" << std::endl
    << "\t--verbose be verbose" << std::endl
    << "\t--clobber clobber the output files" << std::endl
    << "\t--xfactor <f> downsample by a factor in X directon , default 1"<<std::endl
    << "\t--yfactor <f> downsample by a factor in Y directon , default 1"<<std::endl
    << "\t--zfactor <f> downsample by a factor in Z directon , default 2"<<std::endl
    << "\t--factor <f> same as --zfactor"<<std::endl
    << "\t--3dfactor <f>  set X,Y,Z factor to same value"<<std::endl
    << "\t--float store output as float"<<std::endl
    << "\t--short store output as short"<<std::endl
    << "\t--byte store output as byte"<<std::endl;
    

}

int main (int argc, char **argv)
{
  int clobber=0;
  int verbose=0;
  int factor3d=0;
  int store_float=0;
  int store_byte=0;
  int store_short=0;
  
  int factor_x=1;
  int factor_y=1;
  int factor_z=2;
  
  // read the arguments
	static struct option long_options[] =
	  {
		  {"verbose", no_argument, &verbose, 1},
		  {"quiet", no_argument, &verbose, 0},
		  {"clobber", no_argument, &clobber, 1},
      {"xfactor", required_argument, 0, 'x'},
      {"yfactor", required_argument, 0, 'y'},
      {"zfactor", required_argument, 0, 'z'},
      {"factor", required_argument, 0, 'z'},
      {"3dfactor", required_argument, 0, 'd'},
      {"float", no_argument, &store_float, 1},
      {"short", no_argument, &store_short, 1},
      {"byte", no_argument, &store_byte, 1},
		  {0, 0, 0, 0}
	  };

	int c;
  char *_history = time_stamp(argc, argv);
  std::string history=_history;
  free(_history);

	for (;;)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "vf", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
    case 'x':
      factor_x=atoi(optarg);
      if(factor_x<1) 
      {
        std::cerr<<"Error! factor should be >= 1!"<<std::endl;
        return 1;
      }
      break;
    case 'y':
      
      factor_y=atoi(optarg);
      if(factor_y<1) 
      {
        std::cerr<<"Error! factor should be >= 1!"<<std::endl;
        return 1;
      }
      break;
    case 'z':
      factor_z=atoi(optarg);
      if(factor_z<1) 
      {
        std::cerr<<"Error! factor should be >= 1!"<<std::endl;
        return 1;
      }
      break;
		case 'd':
      factor3d=atoi(optarg);
      if(factor3d<2) 
      {
        std::cerr<<"Error! 3d factor should be >= 2!"<<std::endl;
        return 1;
      }
			break;
		case '?':
			/* getopt_long already printed an error message. */
		default:
			show_usage(argv[0]);
			return 1;
		}
	}

	if ((argc - optind) < 2)
	{
		show_usage(argv[0]);
		return 1;
	}
  std::string input=argv[optind];
  std::string output=argv[optind+1];
  
  if (!clobber && !access (output.c_str(), F_OK))
  {
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
  try
  {
    if(factor3d>0)
      factor_x=factor_y=factor_z=factor3d;
    
    minc_1_reader rdr;
    rdr.open(input.c_str());
    
    if((rdr.dim_no()) != 3){
      std::cerr << "Error: volume in " << input.c_str() << " does not have three dimensions." << endl;
      return 1;
    }
    
    simple_volume<float> volume_in;
    load_simple_volume<float>(rdr,volume_in);
    
    simple_volume<float>::idx dims(volume_in.size());
    
    dims/=IDX<size_t>(factor_x,factor_y,factor_z);
    for(int i=0;i<3;i++)
      if(dims[i]<1) dims[i]=1;
    
    simple_volume<float> volume_out(dims);
    
   
    for(int z=0;z<dims[2];z++)
    {
      for(int y=0;y<dims[1];y++)
      for(int x=0;x<dims[0];x++)
      {
        double av=0;
        int c=0;
        
        for(int k=0;k<factor_z;k++)
          for(int j=0;j<factor_y;j++)
            for(int i=0;i<factor_x;i++)
        {
          if(volume_in.hit(x*factor_x+i,
             y*factor_y+j,
             z*factor_z+k))
          {
            av+=volume_in.get(x*factor_x+i,
                                   y*factor_y+j,
                                   z*factor_z+k);
            c++;
          }
        }
        
        if(c==0)
          volume_out.set(x,y,z,0);
        else
          volume_out.set(x,y,z,av/c);
      }
    }
    
    minc::minc_info new_info=rdr.info();
    
    new_info[rdr.map_space(3)].start-=new_info[rdr.map_space(3)].step/2;
    new_info[rdr.map_space(3)].step*=factor_z;
    new_info[rdr.map_space(3)].start+=new_info[rdr.map_space(3)].step/2;    
    new_info[rdr.map_space(3)].length/=factor_z;
    
    new_info[rdr.map_space(2)].start-=new_info[rdr.map_space(2)].step/2;
    new_info[rdr.map_space(2)].step*=factor_y;
    new_info[rdr.map_space(2)].start+=new_info[rdr.map_space(2)].step/2;    
    new_info[rdr.map_space(2)].length/=factor_y;
    
    new_info[rdr.map_space(1)].start-=new_info[rdr.map_space(1)].step/2;
    new_info[rdr.map_space(1)].step*=factor_x;
    new_info[rdr.map_space(1)].start+=new_info[rdr.map_space(1)].step/2;    
    new_info[rdr.map_space(1)].length/=factor_x;
    
    for(int i=1;i<4;i++)
      if(new_info[rdr.map_space(i)].length<1)
        new_info[rdr.map_space(i)].length=1;
    
    minc_1_writer wrt;
    wrt.open(output.c_str(),new_info,2,store_float?NC_FLOAT:store_short?NC_SHORT:store_byte?NC_BYTE:rdr.datatype());
    wrt.copy_headers(rdr);
    wrt.append_history(history.c_str());
    save_simple_volume<float>(wrt,volume_out);
    
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg() << std::endl;
    return 1;
  }
  return 0;
}
