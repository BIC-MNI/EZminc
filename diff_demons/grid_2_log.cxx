#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDisplacementToVelocityFieldLogFilter.h"
#include "itkExponentialDeformationFieldImageFilter.h"
#include "itkDisplacementFieldCompositionFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include <iostream>
#include "mincUtils.h"
#include <getopt.h>

template < class TPixel = float, unsigned int VImageDimension = 3 >

class CommandIterationUpdate: public itk::Command
{
  public:
    typedef  CommandIterationUpdate    Self;
    typedef  itk::Command              Superclass;
    typedef  itk::SmartPointer< Self > Pointer;

    typedef itk::Image< TPixel, VImageDimension >           InternalImageType;
    typedef itk::Vector< TPixel, VImageDimension >          VectorPixelType;
    typedef itk::Image<  VectorPixelType, VImageDimension > VelocityFieldType;
    typedef itk::Image<  VectorPixelType, VImageDimension > DeformationFieldType;
    typedef itk::DisplacementFieldCompositionFilter< DeformationFieldType, DeformationFieldType >
    CompositionFilterType;

    typedef itk::DisplacementToVelocityFieldLogFilter< DeformationFieldType, VelocityFieldType >
    VelocitorFilterType;

    typedef itk::ExponentialDeformationFieldImageFilter< VelocityFieldType, DeformationFieldType >
    ExponentiatorFilterType;

    itkNewMacro ( Self );

    void Execute ( itk::Object *caller, const itk::EventObject & event )
    {
      Execute ( ( const itk::Object * ) caller, event );
    }

    void Execute ( const itk::Object *object, const itk::EventObject & event )
    {
      if ( ! ( itk::IterationEvent().CheckEvent ( &event ) ) )
      {
        return;
      }

      if ( VelocitorFilterType * filter = dynamic_cast< VelocitorFilterType * > (
                                            const_cast< Object * > ( object ) ) )
      {
        typename VelocityFieldType::Pointer velField = filter->GetOutput();

        typename ExponentiatorFilterType::Pointer exponentiator = ExponentiatorFilterType::New();
        exponentiator->SetInput ( velField );

        try
        {
          exponentiator->Update();
        }
        catch ( itk::ExceptionObject & e )
        {
          std::cerr << e;
          itkExceptionMacro ( << "Error in exponentiation" );
        }

        typedef itk::ImageRegionIterator< DeformationFieldType >      IteratorType;

        typedef itk::ImageRegionConstIterator< DeformationFieldType > ConstIteratorType;
        ConstIteratorType itIn ( filter->GetInput(),
                                 filter->GetInput()->GetLargestPossibleRegion() );
        IteratorType      it ( exponentiator->GetOutput(),
                               exponentiator->GetOutput()->GetLargestPossibleRegion() );

        double mse = 0.0;

        double maxError = -1.0;
        double minError = 1e9;

        while ( !it.IsAtEnd() )
        {
          double error = ( it.Value() - itIn.Value() ).GetSquaredNorm();

          if ( error < minError )
          {
            minError = error;
          }

          if ( error > maxError )
          {
            maxError = error;
          }

          mse += error;

          ++it;
          ++itIn;
        }

        std::cout << "Elapsed iterations: " << filter->GetElapsedIterations();

        std::cout << " MSE: " << mse / filter->GetInput()->GetLargestPossibleRegion().GetNumberOfPixels()
                  << " min: " << minError << " max: " << maxError << std::endl;
      }
    }

  protected:
    CommandIterationUpdate()
    {}

  private:
};



void show_usage ( const char *name )
{
	std::cerr
	    << "Grid manipulation program, based on the code from http://hdl.handle.net/10380/3060 " << std::endl
      << "Usage: "<<name<<"  <input> <output>" << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--clobber clobber the output files" << std::endl
      << "\t--iter <n> number of iterations , default 10 " << std::endl
      << "\t--sigma <f> smoothing sigma, default 2 " << std::endl
      << "\t--approx <n> number of BCH approximation terms, default 3" << std::endl
      << "\t--steps <n> number of integration steps , default 500" << std::endl 
      << "\t--byte store in bytes " << std::endl
      << "\t--short store in shorts " << std::endl 
      << "\t--float store in floats " << std::endl
      << "\t--exp - compute exponent (inverse of velocity) " << std::endl
      ;
}

int main ( int argc, char *argv[] )
{
	int clobber=0;
  int store_float=0;
  int store_short=0;
  int store_byte=0;
	int verbose=0;

	int calc_exp=0;
	//parameters
	int iter=10;
	int approx=3;
	int steps=500;
	double sigma=2.0;
	
	
	char *_history = time_stamp(argc, argv); 
	std::string minc_history=_history;
  free ( _history );
	
	static struct option long_options[] =
	{
		{"verbose", no_argument, &verbose, 1},
		{"quiet", no_argument, &verbose, 0},
		{"clobber", no_argument, &clobber, 1},
		{"exp", no_argument, &calc_exp, 1},
		{"sigma", required_argument, 0, 's'},
		{"iter", required_argument, 0, 'i'},
		{"approx", required_argument, 0, 'a'},
		{"steps", required_argument, 0, 't'},
    {"float",   no_argument, &store_float, 1},
    {"short",   no_argument, &store_short, 1},
    {"byte",    no_argument, &store_byte,  1},
    {0, 0, 0, 0}
	};

	int c;
	for ( ;; )
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long ( argc, argv, "s:i:a:t:", long_options, &option_index );

		/* Detect the end of the options. */
		if ( c == -1 )
			break;

		switch ( c )
		{
			case 0:
				break;
			case 's':
				sigma=atof(optarg);
				break;
      case 'i':
        iter=atoi(optarg);
        break;
      case 'a':
        approx=atoi(optarg);
        break;
      case 't':
        steps=atoi(optarg);
        break;
      case '?':
				/* getopt_long already printed an error message. */
			default:
				show_usage ( argv[0] );
				return 1;
		}
	}

	if ( ( argc - optind ) < 2 )
	{
		show_usage ( argv[0] );
		return 1;
	}

	std::string input_f=argv[optind];
  std::string output_f=argv[optind+1];
  //Add Minc support
  itk::RegisterMincIO();
	
  if (!clobber && !access (output_f.c_str (), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }
  
	try
	{
		
  typedef itk::Vector< double, 3 >    VectorType;
  typedef itk::Image< VectorType, 3 > ImageType;

  typedef itk::ImageFileReader< ImageType >                                   ReaderType;
  typedef itk::ImageFileWriter< ImageType >                                   WriterType;
  typedef itk::DisplacementToVelocityFieldLogFilter< ImageType, ImageType >   VelocitorFilterType;
  typedef itk::ExponentialDeformationFieldImageFilter< ImageType, ImageType > ExponentiatorFilterType;
  typedef itk::DisplacementFieldCompositionFilter< ImageType, ImageType >     CompositionFilterType;

  if(verbose) 
		std::cout << "Reading...";
	
  ImageType::Pointer input_field = 0;
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName ( input_f.c_str() );

		reader->Update();

    input_field = reader->GetOutput();

    input_field->DisconnectPipeline();
  }

	if(verbose)
		std::cout << "Done." << std::endl;

	if(verbose)
		std::cout << "Calculating..."<<std::flush;
	
  ImageType::Pointer output_field = 0;
	
	if(calc_exp)
	{
     ExponentiatorFilterType::Pointer exponentiator = ExponentiatorFilterType::New();
     exponentiator->SetInput ( input_field );
     //exponentiator->ComputeInverseOn();
 
     exponentiator->Update();
 
     output_field = exponentiator->GetOutput();
 
     output_field->DisconnectPipeline();
	} else {
    VelocitorFilterType::Pointer velocitor = VelocitorFilterType::New();
    velocitor->SetInput ( input_field );
    velocitor->SetNumberOfIterations ( iter );
		if(sigma>0.0)
		{
			velocitor->SmoothVelocityFieldOn();
			velocitor->SetSigma ( sigma );
		} else {
			velocitor->SmoothVelocityFieldOff();
		}
		
    velocitor->SetNumberOfExponentialIntegrationSteps ( steps );
    velocitor->SetNumberOfBCHApproximationTerms ( approx );

		CommandIterationUpdate< double, 3 >::Pointer observer = CommandIterationUpdate< double, 3 >::New();
		
		if(verbose)
			velocitor->AddObserver ( itk::IterationEvent(), observer );

		velocitor->Update();

    output_field = velocitor->GetOutput();

    output_field->DisconnectPipeline();
  }

	if(verbose)
		std::cout << "Done." << std::endl;

	minc::copy_metadata(output_field,input_field);	
  if(!minc_history.empty())
    minc::append_history ( output_field, minc_history );
	
	if(store_byte)
		minc::set_minc_storage_type ( output_field, NC_BYTE, false );
	if(store_short)
		minc::set_minc_storage_type ( output_field, NC_SHORT, false );
	if(store_float)
		minc::set_minc_storage_type ( output_field, NC_FLOAT, true );
	
	if(verbose)
		std::cout << "Saving...";
  {
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName ( output_f.c_str() );
    writer->SetInput ( output_field );

      writer->Update();
  }

	}
	catch ( itk::ExceptionObject & e )
	{
		std::cerr << e;
		return -1;
	}

  return EXIT_SUCCESS;
}
