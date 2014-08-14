/* ----------------------------- MNI Header -----------------------------------
@NAME       :  volume_2_csv
@DESCRIPTION:  an example of converting a minc file volume to a text format
@COPYRIGHT  :
              Copyright 2011 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <valarray>
#include <math.h>
#include <limits>
#include <unistd.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>


#ifdef HAVE_MINC4ITK
#include <time_stamp.h>    // for creating minc style history entry
#include "itkMincHelpers.h"
#include "itkMincImageIO.h"
#include "itkMincHelpers.h"
#else
#include "itk4MincHelpers.h"
#include "itkMINCTransformAdapter.h"
#endif //HAVE_MINC4ITK


using namespace  std;
using namespace  minc;

void show_usage (const char * prog)
{
  std::cerr<<"Usage:"<<prog<<" in.mnc out.csv --mask <mask> --clobber --terse"<<endl; 
}



int main (int argc, char **argv)
{
  std::string mask_f;
  int verbose=0,clobber=0;
  int terse=0;
  static struct option long_options[] = {
		{"verbose", no_argument,       &verbose, 1},
		{"quiet",   no_argument,       &verbose, 0},
		{"clobber", no_argument,       &clobber, 1},
		{"mask",    required_argument,   0, 'm'},
    {"terse",   no_argument,       &terse, 1},
		{0, 0, 0, 0}
		};
  
  for (;;) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "m:vd:k:i:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1) break;

      switch (c)
			{
			case 0:
				break;
			case 'v':
				cout << "Version: 0.1" << endl;
				return 0;
      case 'm':
        mask_f=optarg; break;
			case '?':
				/* getopt_long already printed an error message. */
			default:
				show_usage (argv[0]);
				return 1;
			}
    }

	if((argc - optind) < 2) {
		show_usage (argv[0]);
		return 1;
	}
  if (!clobber && !access(argv[optind+1], F_OK))
  {
    cerr << argv[optind+1] << " Exists!" << endl;
    return 1;
  }
  
	try
  {
		minc::mask3d::Pointer  mask(minc::mask3d::New());
    
    //creating a minc reader
    
    itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
    //initializing the reader
    reader->SetFileName(argv[optind]);
    reader->Update();
    
    minc::image3d::Pointer img=reader->GetOutput();
    
    
    if(!mask_f.empty())
    {
      itk::ImageFileReader<minc::mask3d >::Pointer reader = itk::ImageFileReader<minc::mask3d >::New();
      //initializing the reader
      reader->SetFileName(mask_f.c_str());
      reader->Update();
      mask=reader->GetOutput();
    }
    const minc::mask3d::RegionType& reg=mask->GetLargestPossibleRegion();
  
    ofstream out(argv[optind+1]);
    if(out.bad())
      ITK_REPORT_ERROR ("can't open file");

    image3d_iterator it(img,img->GetLargestPossibleRegion());
    if(!terse) out<<"x,y,z,I"<<endl;
    if(out.bad())
      ITK_REPORT_ERROR ("can't write to file");
    
    mask3d::IndexType idx;
    
    for(it.GoToBegin();!it.IsAtEnd();++it) {
      tag_point p;
      img->TransformIndexToPhysicalPoint(it.GetIndex(),p);
      
      if(!mask_f.empty())
      {
        mask->TransformPhysicalPointToIndex(p,idx);
        if(!reg.IsInside(idx) || !mask->GetPixel(idx))
            continue;
      }
      
      if(terse)
        out<<it.Value()<<endl;
      else
      out<<p[0]<<","<<p[1]<<","<<p[2]<<","<<it.Value()<<endl;
      if(out.bad())
        ITK_REPORT_ERROR ("can't write to file");
    }
		return 0;
	} 
#ifdef HAVE_MINC4ITK    
  catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  }
#endif
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  }
  return 0;
};
