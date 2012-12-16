#include "minc_wrappers.h"
//#include "minc_io.h"
#include <stdlib.h>
#include <iostream>
#include <getopt.h>
#include <unistd.h>

#include "strtok.h"
#include <time_stamp.h>    // for creating minc style history entry

#include <itkBinaryBallStructuringElement.h>

//Bimodal T threshold 
#include <itkOtsuThresholdImageFilter.h>
//itk morphology filters
#include <itkBinaryDilateImageFilter.h>
#include <itkFastIncrementalBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>

#include <itkBinaryThinningImageFilter.h>
#include <itkBinaryPruningImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>


#include <itkInvertIntensityImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkBinaryMedianImageFilter.h>

using namespace  std;
using namespace  minc;


typedef itk::ImageToImageFilter<
                    minc::mask3d,
                    minc::mask3d > ImageFilter;

typedef itk::BinaryBallStructuringElement<
                    minc::mask3d::PixelType,
                    3  >  BallStructuringElementType;
                    
typedef itk::BinaryDilateImageFilter<
                    minc::mask3d,
                    minc::mask3d,
                    BallStructuringElementType >  DilateFilterType;

typedef itk::BinaryErodeImageFilter<
                    minc::mask3d,
                    minc::mask3d,
                    BallStructuringElementType >  ErodeFilterType;
                    
typedef itk::OtsuThresholdImageFilter< 
                    minc::image3d, 
                    minc::mask3d >            OtsuThresholdFilter;

typedef itk::BinaryThresholdImageFilter< 
                    minc::image3d, 
                    minc::mask3d >            BinaryThresholdFilter;
                    
typedef itk::BinaryPruningImageFilter <
                    minc::mask3d,
                    minc::mask3d>    BinaryPruningFilter;              

typedef itk::BinaryThinningImageFilter<
                    minc::mask3d,
                    minc::mask3d>    BinaryThiningFilter;

typedef itk::ConnectedComponentImageFilter< 
                    minc::mask3d,
                    minc::mask3d>    BinaryConnectedFilter;                  

typedef itk::RelabelComponentImageFilter<
                    minc::mask3d,
                    minc::mask3d>    BinaryRelabelFilter;

typedef itk::BinaryThresholdImageFilter< 
                    minc::mask3d, 
                    minc::mask3d >     BinaryThresholdMaskFilter;

typedef itk::BinaryMedianImageFilter<                  
                    minc::mask3d, 
                    minc::mask3d > BinaryMedianFilter;
                    
typedef itk::InvertIntensityImageFilter< 
                    minc::mask3d, 
                    minc::mask3d > BinaryInvertFilter;
                    
class op_add {
  public:
  minc_mask_voxel _v;

  op_add(minc_mask_voxel v=0):_v(v)
  {
  }
  
  bool operator!=(const op_add& a)
  {
    return a._v!=_v;
  }
  
  bool operator==(const op_add& a)
  {
    return a._v==_v;
  }
  
  minc_mask_voxel operator()(minc_mask_voxel l)
  {
    return _v+l;
  }
};                      
                    
typedef itk::UnaryFunctorImageFilter< 
                    minc::mask3d,
                    minc::mask3d,
                    op_add>     BinaryAddFilter;
                    
void show_usage (const char *name)
{
  std::cerr 
    << "Usage: "<<name<<" <input> <output> " << endl
    << "--verbose be verbose "    << endl
    << "--clobber clobber output files"<<endl
    << "--exp <operations>"<<endl
    << "\tD[<n>] dilate"<<endl
    << "\tE[<n>] erode"<<endl
    << "\tM[<n>] median"<<endl
    << "\tI[0] invert"<<endl
    << "--threshold <threshold>"<<endl
    << "--bimodal"<<endl;
}

ImageFilter::Pointer construct(mask3d::Pointer input, const std::string& par)
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
        flt->SetDilateValue( 1 );
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
        flt->SetErodeValue( 1 );
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      case 'P': //Prune
      {
        int arg=1;
        if(args.empty() || args[0]<=0) arg=1;
        else arg= args[0];
        cout<<"Prune "<<endl;
        BinaryPruningFilter::Pointer flt=BinaryPruningFilter::New();
        //flt->SetErodeValue( 1 );
        flt->SetIteration(arg);
        //TODO:finish this
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      case 'T': //Thinning
      {
        /*if(args.empty() || args[0]<=0) arg=1;
        else arg= args[0];*/
        cout<<"Thinning "<<endl;
        BinaryThiningFilter::Pointer flt=BinaryThiningFilter::New();
        //flt->SetErodeValue( 1 );
        //TODO:finish this
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      
      /*
      case 'C': //connected image filter
      {
        cout<<"Connected "<<endl;
        BinaryConnectedFilter::Pointer flt=BinaryConnectedFilter::New();
        //flt->SetErodeValue( 1 );
        //TODO:finish this
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      case 'R': //connected image filter
      {
        cout<<"Relabeling "<<endl;
        BinaryRelabelFilter::Pointer flt=BinaryRelabelFilter  ::New();
        //flt->SetErodeValue( 1 );
        //TODO:finish this
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }*/
      
      case 'S': //Threshold
      {
        int arg1,arg2;
        if(args.empty() || args[0]<=0) arg1=1;
        else arg1=args[0];
        
        if(args.size()<2 || args[1]<=0) arg2=255;
        else arg2=args[1];
        
        //else arg= args[0];
        cout<<"Threshold "<<arg1<<" "<<arg2<<endl;
        BinaryThresholdMaskFilter::Pointer flt=BinaryThresholdMaskFilter::New();
        flt->	SetLowerThreshold(arg1);
        flt->	SetUpperThreshold(arg2);
        //TODO:finish this
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      
      case 'A': //Add
      {
        int arg1;
        if(args.empty() || args[0]<0) arg1=1;
        else arg1=args[0];
        op_add op(arg1);
        cout<<"Add "<<arg1<<endl;
        BinaryAddFilter::Pointer flt=BinaryAddFilter::New();
        flt->SetFunctor(op);
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
        BinaryMedianFilter::Pointer flt=BinaryMedianFilter::New();
        minc::image3d::SizeType sz;
        sz.Fill(arg1);
        flt->SetRadius(sz);
        flt->SetBackgroundValue(0);
        flt->SetForegroundValue(1);
        flt->SetInput(input);
        return flt.GetPointer();
        break;
      }
      case 'I':
      {
        cout<<"Invert "<<endl;
        BinaryInvertFilter::Pointer flt=BinaryInvertFilter::New();
        flt->SetMaximum(1);
        flt->SetInput(input);
        //std::cout<<flt<<std::endl;
        return flt.GetPointer();
        break;
      }
      default:
        std::cerr<<"Can't decode:"<<par.c_str()<<std::endl;
    }
  }
  return ImageFilter::Pointer();
}

int main (int argc, char **argv)
{
  int verbose=1;
  int clobber=0;
  double threshold=1e10;
  int bimodal=0;
  std::string operations;
  char *history = time_stamp(argc, argv); 
  
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
    itk::RegisterMincIO();

    mask3d::Pointer mask_img( mask3d::New() ),
                    output_img;
    if(bimodal)
    {
      itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
      reader->SetFileName(input.c_str());
      reader->Update();
      minc::image3d::Pointer img=reader->GetOutput();
      
      OtsuThresholdFilter::Pointer flt=OtsuThresholdFilter::New();
      //flt->SetLowerThreshold(threshold);
      flt->SetOutsideValue( 1 );
      flt->SetInsideValue( 0 );
      flt->SetInput(img);
      flt->Update();
      if(verbose) 
        cout<<"BimodalT threshold:"<<flt->GetThreshold()<<endl;
      mask_img=flt->GetOutput();
    } else if(threshold!=1e10) {
      itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
      reader->SetFileName(input.c_str());
      reader->Update();
      minc::image3d::Pointer img=reader->GetOutput();
      
      BinaryThresholdFilter::Pointer flt=BinaryThresholdFilter::New();
      flt->SetLowerThreshold(threshold);
      flt->SetInput(img);
      flt->SetInsideValue(1);
      flt->SetOutsideValue(0);
      flt->Update();
      mask_img=flt->GetOutput();
    } else {
      itk::ImageFileReader<minc::mask3d >::Pointer reader = itk::ImageFileReader<minc::mask3d >::New();
      reader->SetFileName(input.c_str());
      reader->Update();
      mask_img=reader->GetOutput();
      normalize_mask(mask_img);
    }
    //parse through operations
    if(operations.empty())
    {
      output_img=mask_img;
    } else {
      std::vector<ImageFilter::Pointer> filter_list;
      ImageFilter::Pointer last;
      ImageFilter::Pointer flt;
      
      //pcrepp::Pcre delim("[\\s]+");
      //std::vector<std::string> splitted=delim.split(operations);
      std::vector<std::string> splitted;
      stringtok<std::vector<std::string> >(splitted,operations," ");
      
      for(std::vector<std::string>::iterator  A = splitted.begin(); A != splitted.end(); ++A)
      {
        if(last.IsNull())
          flt=construct( mask_img, *A);
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
        output_img=mask_img;
        cout<<"No processing!"<<endl;
      } else {
        last->Update();
        output_img=last->GetOutput();
      }
    }
    minc::copy_metadata(output_img,mask_img);
    minc::append_history(output_img,history);
    free(history);
    
    itk::ImageFileWriter< minc::mask3d >::Pointer writer = itk::ImageFileWriter<minc::mask3d >::New();
    writer->SetFileName(output.c_str());
    
    writer->SetInput( output_img );
    writer->UseInputMetaDataDictionaryOn();
    
    writer->Update();

    
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
  
}
