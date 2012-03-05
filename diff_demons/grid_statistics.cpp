#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDisplacementToVelocityFieldLogFilter.h"
#include <itkExponentialDeformationFieldImageFilter.h>
#include "itkDisplacementFieldCompositionFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <itkMinimumMaximumImageCalculator.h>
#include <itkWarpHarmonicEnergyCalculator.h>
#include <itkDisplacementFieldJacobianDeterminantFilter.h>

#include <iostream>
#include "mincUtils.h"
#include <getopt.h>

void show_usage ( const char *name )
{
  std::cerr
      << "Grid stats program, based on the code from http://hdl.handle.net/10380/3060 " << std::endl
      << "Usage: "<<name<<"  <input field> " << std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--exp  - input field is in log space" << std::endl
      << " Output:"<<std::endl
      << "minJac,Q002,Q01,Q99,Q998,maxJac,harmonic,JacBelow0"<<std::endl;
     ;

}



int main ( int argc, char *argv[] )
{
  int verbose=0;
  int exp_transform=0;
  
  
  std::string mask_f;
  char *_history = time_stamp(argc, argv); 
  std::string minc_history=_history;
  free ( _history );
  
  static struct option long_options[] =
  {
    {"verbose", no_argument, &verbose, 1},
    {"quiet", no_argument, &verbose, 0},
    {"exp", no_argument, &exp_transform, 1},
    {"mask", required_argument, 0, 'm'},
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
      case 'm':
        mask_f=optarg;
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage ( argv[0] );
        return 1;
    }
  }

  if ( ( argc - optind ) < 1 )
  {
    show_usage ( argv[0] );
    return 1;
  }

  std::string input_f=argv[optind];
  //Add Minc support
  
  
  try
  {
    itk::ObjectFactoryBase::RegisterFactory ( itk::MincImageIOFactory::New() );
    
    typedef float TPixel;
    const unsigned int VImageDimension = 3;//TODO: make this a template?
    
    typedef itk::Image< TPixel, VImageDimension >           InternalImageType;
    typedef itk::Vector< TPixel, VImageDimension >          VectorPixelType;
    typedef itk::Image<  VectorPixelType, VImageDimension > VectorImageType;
    typedef itk::Image<  VectorPixelType, VImageDimension > VectorImageType;
    
    typedef itk::DisplacementFieldJacobianDeterminantFilter< VectorImageType, TPixel > JacobianFilterType;
    typedef itk::MinimumMaximumImageCalculator < InternalImageType >     MinMaxFilterType;
    typedef itk::WarpHarmonicEnergyCalculator < VectorImageType >   HarmonicEnergyCalculatorType;
    typedef itk::VectorCentralDifferenceImageFunction < VectorImageType > WarpGradientCalculatorType;
    typedef WarpGradientCalculatorType::OutputType WarpGradientType;
    typedef itk::ExponentialDeformationFieldImageFilter< VectorImageType, VectorImageType > ExponentiatorFilterType;
    typedef itk::ImageFileReader< VectorImageType >  VectorReaderType;


    VectorImageType::Pointer input_vector_field;
    VectorImageType::Pointer input_def_field;;
    
    
    VectorReaderType::Pointer reader = VectorReaderType::New();
    reader->SetFileName ( input_f.c_str() );

    reader->Update();
    input_vector_field = reader->GetOutput();
    input_vector_field->DisconnectPipeline();
    

    if(exp_transform)
    {
      ExponentiatorFilterType::Pointer exponentiator = ExponentiatorFilterType::New();
      exponentiator->SetInput(input_vector_field);
      
      exponentiator->Update();
      
      input_def_field = exponentiator->GetOutput();
      input_def_field->DisconnectPipeline();
    } else {
      input_def_field = input_vector_field; //use it as it is 
    }
    
    
    JacobianFilterType::Pointer m_JacobianFilter = JacobianFilterType::New();
    m_JacobianFilter->SetUseImageSpacing ( true );
    m_JacobianFilter->ReleaseDataFlagOn();

    MinMaxFilterType::Pointer m_Minmaxfilter = MinMaxFilterType::New();

    HarmonicEnergyCalculatorType::Pointer m_HarmonicEnergyCalculator = HarmonicEnergyCalculatorType::New();    
  
    m_HarmonicEnergyCalculator->SetImage ( input_def_field );

    m_HarmonicEnergyCalculator->Compute();
    const double harmonicEnergy= m_HarmonicEnergyCalculator->GetHarmonicEnergy();


    m_JacobianFilter->SetInput ( input_def_field );
    m_JacobianFilter->UpdateLargestPossibleRegion();


    const unsigned int numPix = m_JacobianFilter->
                                GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();

    TPixel* pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
    TPixel* pix_end = pix_start + numPix;

    TPixel* jac_ptr;

    // Get percentage of det(Jac) below 0
    unsigned int jacBelowZero ( 0u );

    for ( jac_ptr = pix_start; jac_ptr != pix_end; ++jac_ptr )
    {
      if ( *jac_ptr <= 0.0 ) ++jacBelowZero;
    }

    const double jacBelowZeroPrc = static_cast<double> ( jacBelowZero ) / static_cast<double> ( numPix );


    // Get min an max jac
    const double minJac = * ( std::min_element ( pix_start, pix_end ) );
    const double maxJac = * ( std::max_element ( pix_start, pix_end ) );

    // Get some quantiles
    // We don't need the jacobian image
    // we can modify/sort it in place
    jac_ptr = pix_start + static_cast<unsigned int> ( 0.002 * numPix );

    std::nth_element ( pix_start, jac_ptr, pix_end );

    const double Q002 = *jac_ptr;

    jac_ptr = pix_start + static_cast<unsigned int> ( 0.01 * numPix );

    std::nth_element ( pix_start, jac_ptr, pix_end );

    const double Q01 = *jac_ptr;

    jac_ptr = pix_start + static_cast<unsigned int> ( 0.99 * numPix );

    std::nth_element ( pix_start, jac_ptr, pix_end );

    const double Q99 = *jac_ptr;

    jac_ptr = pix_start + static_cast<unsigned int> ( 0.998 * numPix );

    std::nth_element ( pix_start, jac_ptr, pix_end );

    const double Q998 = *jac_ptr;
    
    std::cout<<minJac<<","
             <<Q002<<","
             <<Q01<<","
             <<Q99<<","
             <<Q998<<","
             <<Q99<<","
             <<harmonicEnergy<<","
             <<jacBelowZeroPrc<<std::endl;
  }
  catch ( itk::ExceptionObject & e )
  {
    std::cerr << e;
    return -1;
  }

  return EXIT_SUCCESS;
}
