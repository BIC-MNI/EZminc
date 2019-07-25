#include "minc_wrappers.h"
#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <unistd.h>
#include <time_stamp.h>    // for creating minc style history entry


using namespace  std;
using namespace  minc;
                    
void show_usage (const char *name)
{
  std::cerr 
    << "Usage: "<<name<<" <input> <output> " << endl
    << "--verbose be verbose "    << endl
    << "--clobber clobber output files"<<endl;

}

int main (int argc, char **argv)
{
  int verbose=1;
  int clobber=0;
  int inverse=0;
  double threshold=1e10;
  int bimodal=0;
  std::string operations;
  char *history = time_stamp(argc, argv); 
  
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
    {"inverse", no_argument, &inverse, 1},
    {0, 0, 0, 0}
    };
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int
    option_index = 0;

    c = getopt_long (argc, argv, "", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
      break;
    case 'v':
      cout << "Version: 1.0" << endl;
      return 0;
    case '?':
      /* getopt_long already printed an error message. */
    default:
      show_usage (argv[0]);
      return 1;
    }
  }
  
  if((argc - optind)<2)
  {
    show_usage(argv[0]);
    return 1;
  }
  std::string input=argv[optind];
  std::string output=argv[optind+1];
  if (!clobber && !access (output.c_str (), F_OK))
  {
    cerr << output.c_str () << " Exists!" << endl;
    return 1;
  }
  
  try
  {
#ifdef HAVE_MINC4ITK
    itk::RegisterMincIO();
#endif
    
    itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
    reader->SetFileName(input.c_str());
    reader->Update();
    minc::image3d::Pointer img=reader->GetOutput();
    minc::image3d::DirectionType dir_cos=img->GetDirection();
    minc::image3d::PointType origin=img->GetOrigin();
    
    itk::Vector<double,3> translation;
    translation.Fill(0);
    
    if(inverse)
    {
      dir_cos=dir_cos.GetInverse();
    }
    write_linear_xfm(output.c_str(),dir_cos,translation);
    //now write out xfm file
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  }
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
  return 0;
}
