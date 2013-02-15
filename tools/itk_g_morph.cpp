#include <stdlib.h>
#include <iostream>
#include <getopt.h>
#include <unistd.h>

#include "strtok.h"

#include <itkImage.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkMedianImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleErodeImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>

#if ( ITK_VERSION_MAJOR < 4 )
#include <time_stamp.h>    // for creating minc style history entry
#include "itkMincImageIOFactory.h"
#include "itkMincImageIO.h"
#include "itkMincHelpers.h"
using namespace  minc;
#endif

using namespace  std;


typedef itk::Image<unsigned char,3> MaskType;
typedef itk::Image<float,3> ImageType;

typedef itk::ImageToImageFilter<
                   ImageType,
                   ImageType > ImageFilter;

typedef itk::BinaryBallStructuringElement<
                   MaskType::PixelType,
                    3  >  BallStructuringElementType;

typedef itk::MedianImageFilter<
              ImageType,ImageType >  MedianFilterType;


typedef itk::GrayscaleErodeImageFilter<
                           ImageType,
                           ImageType,
                            BallStructuringElementType >  ErodeFilterType;

typedef itk::GrayscaleDilateImageFilter<
                           ImageType,
                           ImageType,
                            BallStructuringElementType >  DilateFilterType;

                    

void show_usage (const char *name)
{
  std::cerr 
	  << "Usage: "<<name<<" <input> <output> " << endl
    << "--verbose be verbose "    << endl
    << "--clobber clobber output files"<<endl
    << "--exp <operations>"<<endl
    << "\tD[<n>] dilate"<<endl
    << "\tE[<n>] erode"<<endl
    << "\tM[<n>] median"<<endl;
}

ImageFilter::Pointer construct(ImageType::Pointer input, const std::string& par)
{
  if(par.empty()) return ImageFilter::Pointer();
  //try 
  {
    //pcrepp::Pcre analyze("(\\S)\\[(.*)\\]");
    //if(!analyze.search(par)) return ImageFilter::Pointer();
    size_t left_bk =par.find('[');
    size_t right_bk=par.find(']');
    
    if(left_bk==string::npos || right_bk==string::npos)
    {
      std::cerr<<"WARNING: Unknown operation:"<<par.c_str()<<std::endl;
      return ImageFilter::Pointer();
    }
    std::string op=par.substr(0,1);
    std::vector<int> args;
    //if(analyze.matches()>1)
    {
      //pcrepp::Pcre split("\\s*,\\s*");
      std::vector<std::string> splitted;
      stringtok<std::vector<std::string> >(splitted,par.substr(left_bk+1,right_bk-left_bk-1));
      for(std::vector<std::string>::iterator A = splitted.begin(); A != splitted.end(); ++A)
      {
        args.push_back(atoi((*A).c_str()));
      }
    }
    
    switch(op[0])
    {
      case 'D': //Dilate
      {
        int arg=1;
        if(!args.empty() && args[0]>0) arg= args[0];
        cout<<"Dilate "<<arg<<endl;
        DilateFilterType::Pointer flt=DilateFilterType::New();
        BallStructuringElementType  structuringElement;
        //unsigned long rad[3]={ arg, arg, arg};
        structuringElement.SetRadius( arg );
        structuringElement.CreateStructuringElement();
        flt->SetKernel( structuringElement );
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      
      case 'E': //Erode
      {
        int arg=1;
        if(args.empty() || args[0]<=0) arg=1;
        else arg= args[0];
        cout<<"Erode "<<arg<<endl;
        ErodeFilterType::Pointer flt=ErodeFilterType::New();
        BallStructuringElementType  structuringElement;

        structuringElement.SetRadius( arg );
        structuringElement.CreateStructuringElement();
        flt->SetKernel( structuringElement );
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      
      case 'M': //Median
      {
        int arg1;
        if(args.empty() || args[0]<0) arg1=1;
        else arg1=args[0];
        cout<<"Median "<<arg1<<endl;
        MedianFilterType::Pointer flt=MedianFilterType::New();
       ImageType::SizeType sz;
        sz.Fill(arg1);
        flt->SetRadius(sz);
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      default:
        std::cerr<<"Can't decode:"<<par.c_str()<<std::endl;
    }
  }
  //catch (pcrepp::Pcre::exception &E) {
   /*
    * the Pcre class has thrown an exception
    */
   //cerr << "Pcre++ error: " << E.what() << endl;
  //}
  
  return ImageFilter::Pointer();
}

int main (int argc, char **argv)
{
  int verbose=1;
  int clobber=0;
  double threshold=1e10;
  int bimodal=0;
  std::string operations;
#if ( ITK_VERSION_MAJOR < 4 )
  char *history = time_stamp(argc, argv); 
#endif
  
  static struct option long_options[] = { 
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"clobber", no_argument, &clobber, 1},
		{"exp",  required_argument, 0, 'e'},
    {"threshold",  required_argument, 0, 't'},
    {"bimodal", no_argument, &bimodal, 1},
		{0, 0, 0, 0}
		};
	int c;
	for (;;)
	{
		/* getopt_long stores the option index here. */
		int
		option_index = 0;

		c = getopt_long (argc, argv, "s:t:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			break;
		case 'e':
			operations = optarg;
			break;
		case 't':
			threshold = atof(optarg);
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
#if ( ITK_VERSION_MAJOR < 4 )
    itk::RegisterMincIO();
#endif
    
    itk::ImageFileReader<ImageType >::Pointer reader = itk::ImageFileReader<ImageType >::New();
    
    ImageType::Pointer output_img;
    
    reader->SetFileName(input.c_str());
    reader->Update();
    
    ImageType::Pointer img=reader->GetOutput();
    
    //parse through operations
    if(operations.empty())
    {
      output_img=img;
    } else {
      std::vector<ImageFilter::Pointer> filter_list;
      ImageFilter::Pointer last;
      ImageFilter::Pointer flt;
      
      std::vector<std::string> splitted;
      stringtok<std::vector<std::string> >(splitted,operations," ");
      
      for(std::vector<std::string>::iterator  A = splitted.begin(); A != splitted.end(); ++A)
      {
        if(last.IsNull())
          flt=construct( img, *A);
        else
          flt=construct( last->GetOutput(), *A);
        
        if(!flt.IsNull())
        {
          filter_list.push_back(flt);
          last=flt;
        }
        else
          cout<<"Zero output!"<<endl;
      }
      
      if(last.IsNull())
      {
        output_img=img;
        cout<<"No processing!"<<endl;
      } else {
        last->Update();
        output_img=last->GetOutput();
      }
    }
#if ( ITK_VERSION_MAJOR < 4 )
    minc::copy_metadata(output_img,img);
    minc::append_history(output_img,history);
    free(history);
#endif 
    
    itk::ImageFileWriter< ImageType >::Pointer writer = itk::ImageFileWriter<ImageType >::New();
    writer->SetFileName(output.c_str());
    
    writer->SetInput( output_img );
    writer->UseInputMetaDataDictionaryOn();
    
    writer->Update();
    
  } 
#if ( ITK_VERSION_MAJOR < 4 )
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
  
}
