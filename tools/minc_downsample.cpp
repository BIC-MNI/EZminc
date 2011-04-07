#include <iostream>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>

#include <time_stamp.h>    // for creating minc style history entry
// for get_opt_long
#include <getopt.h>
#include <stdlib.h>

using namespace  std;
using namespace  minc;


void show_usage(const char *name)
{
  std::cerr 
	  << "Usage: "<<name<<" <input> <output> " << std::endl
    << "Optional parameters:" << std::endl
    << "\t--verbose be verbose" << std::endl
    << "\t--clobber clobber the output files" << std::endl
    << "\t--factor <f> downsample by a factor , default 2"<<std::endl
    << "\t--3dfactor <f> downsample by a 3dfactor , default 1"<<std::endl
    << "\t--float store output as float"<<std::endl
    << "\t--short store output as short"<<std::endl
    << "\t--byte store output as byte"<<std::endl;
    

}

int main (int argc, char **argv)
{
  int clobber=0;
  int verbose=0;
  int factor=2;
  int factor3d=1;
  int store_float=0;
  int store_byte=0;
  int store_short=0;
  // read the arguments
	static struct option long_options[] =
	  {
		  {"verbose", no_argument, &verbose, 1},
		  {"quiet", no_argument, &verbose, 0},
		  {"clobber", no_argument, &clobber, 1},
      {"factor", required_argument, 0, 'f'},
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
		case 'f':
      factor=atoi(optarg);
      if(factor<2) 
      {
        std::cerr<<"Error! factor should be >= 2!"<<std::endl;
        return 1;
      }
			break;
		case 'd':
      factor3d=factor=atoi(optarg);
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
    int factor_x=factor3d;
    int factor_y=factor3d;
    int factor_z=factor;
    
    minc_1_reader rdr;
    rdr.open(input.c_str());
    
    if((rdr.dim_no()) != 3){
      std::cerr << "Error: volume in " << input.c_str() << " does not have three dimensions." << endl;
      return 1;
    }
    
    simple_volume<float> volume_in;
    load_simple_volume<float>(rdr,volume_in);
    
    simple_volume<float> volume_out(volume_in.dim(0)/factor_x,
                                    volume_in.dim(1)/factor_y,
                                    volume_in.dim(2)/factor_z);
    
   
    for(int z=0;z<volume_out.dim(2);z++)
    {
      for(int y=0;y<volume_out.dim(1);y++)
      for(int x=0;x<volume_out.dim(0);x++)
      {
        double av=0;
        int c=0;
        
        for(int k=0;k<factor_z;k++)
          for(int j=0;j<factor_y;j++)
            for(int i=0;i<factor_x;i++)
        {
            av+=volume_in.safe_get(x*factor_x+i,
                                   y*factor_y+j,
                                   z*factor_z+k);
            c++;
        }
        
        if(!c) continue;
        
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
