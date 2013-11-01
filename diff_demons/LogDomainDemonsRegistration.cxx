/**
 * This is a small tool that shows how to use the log-domain demons algorithm.
 * The user can choose if standard or symmetrized log-domain demons should be used.
 * The user can also choose the type of demons forces, or other parameters.
 *
 * \author Florence Dru, INRIA and Tom Vercauteren, MKT
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

#include <itkCommand.h>
#include <itkLogDomainDemonsRegistrationFilter.h>
#include <itkSymmetricLogDomainDemonsRegistrationFilter.h>
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <itkGridForwardWarpImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiResolutionLogDomainDeformableRegistration.h>
#include <itkTransformFileReader.h>
#include <itkTransformToVelocityFieldSource.h>
#include <itkVectorCentralDifferenceImageFunction.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkWarpHarmonicEnergyCalculator.h>
#include <itkWarpImageFilter.h>


#if ( ITK_VERSION_MAJOR > 3 ) 
#include <itkTransformToDisplacementFieldSource.h>
#else
#include <itkTransformToDeformationFieldSource.h>
#endif

#include <itkDisplacementToVelocityFieldLogFilter.h>

#include <metaCommand.h>

#include <errno.h>
#include <iostream>
#include <limits.h>

//MINC support
#include "mincUtils.h"

extern std::string minc_history;

struct arguments
{
  std::string  fixedImageFile;  /* -f option */
  std::string  movingImageFile; /* -m option */
  std::string  fixedImageMaskFile;  /* --fixed-mask option */
  std::string  movingImageMaskFile; /* --moving-mask option */
  std::string  inputFieldFile;  /* -b option */
  std::string  inputTransformFile;  /* -p option */
  std::string  outputImageFile; /* -o option */
  std::string  outputJacobianFile; /* -J option */
  std::string  outputDeformationFieldFile; /* -O */
  std::string  outputInverseDeformationFieldFile; /* -I */
  std::string  outputVelocityFieldFile; /* -V */
  std::string  trueFieldFile;   /* -r option */
  std::vector<unsigned int> numIterations;   /* -i option */
  float sigmaVel;               /* -s option */
  float sigmaUp;                /* -g option */
  float maxStepLength;          /* -l option */
  unsigned int updateRule;      /* -a option */
  unsigned int gradientType;    /* -t option */
  unsigned int NumberOfBCHApproximationTerms; /* -c option */
  bool useHistogramMatching;    /* -e option */
  unsigned int verbosity;       /* -d option */

  friend std::ostream& operator<< ( std::ostream& o, const arguments& args )
  {
    std::ostringstream osstr;

    for ( unsigned int i = 0; i < args.numIterations.size(); ++i )
    {
      osstr << args.numIterations[i] << " ";
    }

    std::string iterstr = "[ " + osstr.str() + "]";

    std::string gtypeStr;

    switch ( args.gradientType )
    {
      case 0:
        gtypeStr = "symmetrized (ESM for diffeomorphic and compositive)";
        break;
      case 1:
        gtypeStr = "fixed image (Thirion's vanilla forces)";
        break;
      case 2:
        gtypeStr = "warped moving image (Gauss-Newton for diffeomorphic and compositive)";
        break;
      case 3:
        gtypeStr = "mapped moving image (Gauss-Newton for additive)";
        break;
      default:
        gtypeStr = "unsuported";
    }

    std::string uruleStr;

    switch ( args.updateRule )
    {
      case 0:
        uruleStr = "BCH approximation on velocity fields (log-domain)";
        break;
      case 1:
        uruleStr = "Symmetrized BCH approximation on velocity fields (symmetric log-domain)";
        break;
      default:
        uruleStr = "unsuported";
    }

    std::string histoMatchStr = ( args.useHistogramMatching ? "true" : "false" );

    return o
           << "Arguments structure:" << std::endl
           << "  Fixed image file: " << args.fixedImageFile << std::endl
           << "  Moving image file: " << args.movingImageFile << std::endl
           << "  Fixed image mask file: " << args.fixedImageMaskFile << std::endl
           << "  Moving image mask file: " << args.movingImageMaskFile << std::endl
           << "  Input velocity field file: " << args.inputFieldFile << std::endl
           << "  Input transform file: " << args.inputTransformFile << std::endl
           << "  Output image file: " << args.outputImageFile << std::endl
           << "  Output deformation field file: " << args.outputDeformationFieldFile << std::endl
           << "  Output inverse deformation field file: " << args.outputInverseDeformationFieldFile << std::endl
           << "  Output velocity field file: " << args.outputVelocityFieldFile << std::endl
           << "  Output jacobian field file: " << args.outputJacobianFile << std::endl
           << "  True deformation field file: " << args.trueFieldFile << std::endl
           << "  Number of multiresolution levels: " << args.numIterations.size() << std::endl
           << "  Number of log-domain demons iterations: " << iterstr << std::endl
           << "  Velocity field sigma: " << args.sigmaVel << std::endl
           << "  Update field sigma: " << args.sigmaUp << std::endl
           << "  Maximum update step length: " << args.maxStepLength << std::endl
           << "  Update rule: " << uruleStr << std::endl
           << "  Type of gradient: " << gtypeStr << std::endl
           << "  Number of terms in the BCH expansion: " << args.NumberOfBCHApproximationTerms << std::endl
           << "  Use histogram matching: " << histoMatchStr << std::endl
           << "  Algorithm verbosity (debug level): " << args.verbosity;
  }
};

void help_callback()
{
  std::cout << std::endl;
	std::cout << "Code based on http://hdl.handle.net/10380/3060 " << std::endl;
  std::cout << "Copyright (c) 2009 INRIA and Mauna Kea Technologies" << std::endl;
  std::cout << "Code: Florence Dru and Tom Vercauteren" << std::endl;
  std::cout << "With additions by Vladimir S. FONOV" << std::endl;
  std::cout << "Report bugs to <tom.vercauteren \\at maunakeatech.com>" << std::endl;

  exit ( EXIT_FAILURE );
};

int atoi_check ( const char * str )
{
  char *endptr;
  long val = strtol ( str, &endptr, 0 );

  /* Check for various possible errors */

  if ( ( errno == ERANGE && ( val == LONG_MAX || val == LONG_MIN ) )
       || val >= INT_MAX || val <= INT_MIN )
  {
    std::cout << std::endl;
    std::cout << "Cannot parse integer. Out of bound." << std::endl;
    exit ( EXIT_FAILURE );
  }

  if ( endptr == str || *endptr != '\0' )
  {
    std::cout << std::endl;
    std::cout << "Cannot parse integer. Contains non-digits or is empty." << std::endl;
    exit ( EXIT_FAILURE );
  }

  return val;
}


std::vector<unsigned int> parseUIntVector ( const std::string & str )
{
  std::vector<unsigned int> vect;

  std::string::size_type crosspos = str.find ( 'x', 0 );

  if ( crosspos == std::string::npos )
  {
    // only one uint
    vect.push_back ( static_cast<unsigned int> ( atoi_check ( str.c_str() ) ) );
    return vect;
  }

  // first uint
  vect.push_back ( static_cast<unsigned int> (
                     atoi_check ( ( str.substr ( 0, crosspos ) ).c_str() ) ) );

  while ( true )
  {
    std::string::size_type crossposfrom = crosspos;
    crosspos =  str.find ( 'x', crossposfrom + 1 );

    if ( crosspos == std::string::npos )
    {
      vect.push_back ( static_cast<unsigned int> (
                         atoi_check ( ( str.substr ( crossposfrom + 1, str.length() - crossposfrom - 1 ) ).c_str() ) ) );
      return vect;
    }

    vect.push_back ( static_cast<unsigned int> (

                       atoi_check ( ( str.substr ( crossposfrom + 1, crosspos - crossposfrom - 1 ) ).c_str() ) ) );
  }
}


void parseOpts ( int argc, char **argv, struct arguments & args )
{
  // Command line parser
  MetaCommand command;
  command.SetParseFailureOnUnrecognizedOption ( true );
  command.SetHelpCallBack ( help_callback );

  // Fill some information about the software
  command.SetAuthor ( "Florence Dru and Tom Vercauteren with additions from Vladimir S. FONOV" );

  command.SetAcknowledgments ( "This work is supported by INRIA (Asclepios team) and Mauna Kea Technologies" );

  command.SetDescription ( "Basic image registration tool with the log-domain demons algorithm." );

  // Define parsing options
  command.SetOption ( "FixedImageFile", "f", true, "Fixed image filename" );
  command.SetOptionLongTag ( "FixedImageFile", "fixed-image" );
  command.AddOptionField ( "FixedImageFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "MovingImageFile", "m", true, "Moving image filename" );
  command.SetOptionLongTag ( "MovingImageFile", "moving-image" );
  command.AddOptionField ( "MovingImageFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "FixedImageMaskFile", "", false, "Fixed mask filename" );
  command.SetOptionLongTag ( "FixedImageMaskFile", "fixed-mask" );
  command.AddOptionField ( "FixedImageMaskFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "MovingImageMaskFile", "", false, "Moving mask filename" );
  command.SetOptionLongTag ( "MovingImageMaskFile", "moving-mask" );
  command.AddOptionField ( "MovingImageMaskFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "InputFieldFile", "b", false, "Input velocity field filename" );
  command.SetOptionLongTag ( "InputFieldFile", "input-field" );
  command.AddOptionField ( "InputFieldFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "InputTransformFile", "p", false, "Input transform filename" );
  command.SetOptionLongTag ( "InputTransformFile", "input-transform" );
  command.AddOptionField ( "InputTransformFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "OutputImageFile", "o", false, "Output image filename" );
  command.SetOptionLongTag ( "OutputImageFile", "output-image" );
  command.AddOptionField ( "OutputImageFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "OutputDeformationFieldFile", "O", false, "Output deformation field filename" );
  command.SetOptionLongTag ( "OutputDeformationFieldFile", "outputDef-field" );
  command.AddOptionField ( "OutputDeformationFieldFile", "filename", MetaCommand::STRING, false );

  command.SetOption ( "OutputInverseDeformationFieldFile", "I", false, "Output inverse deformation field filename" );
  command.SetOptionLongTag ( "OutputInverseDeformationFieldFile", "outputInvDef-field" );
  command.AddOptionField ( "OutputInverseDeformationFieldFile", "filename", MetaCommand::STRING, false );

  command.SetOption ( "OutputVelocityFieldFile", "V", false, "Output velocity field filename" );
  command.SetOptionLongTag ( "OutputVelocityFieldFile", "outputVel-field" );
  command.AddOptionField ( "OutputVelocityFieldFile", "filename", MetaCommand::STRING, false );

  command.SetOption ( "OutputJacobianFile", "J", false, "Generate the jacobian of  output field " );
  command.SetOptionLongTag ( "OutputJacobianFile", "output-jacobian" );
  command.AddOptionField ( "OutputJacobianFile", "filename", MetaCommand::STRING, false );
  
  command.SetOption ( "TrueFieldFile", "r", false, "Specify a \"true\" deformation field to compare the registration result with (useful for synthetic experiments)" );
  command.SetOptionLongTag ( "TrueFieldFile", "true-field" );
  command.AddOptionField ( "TrueFieldFile", "filename", MetaCommand::STRING, true );

  command.SetOption ( "NumberOfIterationsPerLevels", "i", false, "List of number of iterations for each multi-scale pyramid level < UINTx...xUINT >" );
  command.SetOptionLongTag ( "NumberOfIterationsPerLevels", "num-iterations" );
  command.AddOptionField ( "NumberOfIterationsPerLevels", "uintvect", MetaCommand::STRING, true, "15x10x5" );

  command.SetOption ( "VelocityFieldSigma", "s", false, "Smoothing sigma for the velocity field (pixel units). Setting it below 0.5 means no smoothing will be performed" );
  command.SetOptionLongTag ( "VelocityFieldSigma", "vel-field-sigma" );
  command.AddOptionField ( "VelocityFieldSigma", "floatval", MetaCommand::FLOAT, true, "1.5" );

  command.SetOption ( "UpdateFieldSigma", "g", false, "Smoothing sigma for the update field (pixel units). Setting it below 0.5 means no smoothing will be performed" );
  command.SetOptionLongTag ( "UpdateFieldSigma", "up-field-sigma" );
  command.AddOptionField ( "UpdateFieldSigma", "floatval", MetaCommand::FLOAT, true, "0.0" );

  command.SetOption ( "MaximumUpdateStepLength", "l", false, "Maximum length of an update vector (pixel units). Setting it to 0 implies no restrictions will be made on the step length" );
  command.SetOptionLongTag ( "MaximumUpdateStepLength", "max-step-length" );
  command.AddOptionField ( "MaximumUpdateStepLength", "floatval", MetaCommand::FLOAT, true, "2.0" );

  command.SetOption ( "UpdateRule", "a", false, "Type of update rule. 0: exp(v) <- exp(v) o exp(u) (log-domain), 1: exp(v) <- symmetrized( exp(v) o exp(u) ) (symmetric log-domain)" );
  command.SetOptionLongTag ( "UpdateRule", "update-rule" );
  command.AddOptionField ( "UpdateRule", "type", MetaCommand::INT, true, "1" );
  command.SetOptionRange ( "UpdateRule", "type", "0", "1" );

  command.SetOption ( "GradientType", "t", false, "Type of gradient used for computing the demons force. 0 is symmetrized, 1 is fixed image, 2 is warped moving image, 3 is mapped moving image" );
  command.SetOptionLongTag ( "GradientType", "gradient-type" );
  command.AddOptionField ( "GradientType", "type", MetaCommand::INT, true, "0" );
  command.SetOptionRange ( "GradientType", "type", "0", "3" );

  command.SetOption ( "NumberOfBCHApproximationTerms", "c", false, "Number of terms in the BCH expansion" );
  command.SetOptionLongTag ( "NumberOfBCHApproximationTerms", "num-bch-terms" );
  command.AddOptionField ( "NumberOfBCHApproximationTerms", "intval", MetaCommand::INT, true, "2" );
  command.SetOptionRange ( "NumberOfBCHApproximationTerms", "intval", "2", "4" );

  command.SetOption ( "UseHistogramMatching", "e", false, "Use histogram matching prior to registration (e.g. for different MR scanners)" );
  command.SetOptionLongTag ( "UseHistogramMatching", "use-histogram-matching" );
  command.AddOptionField ( "UseHistogramMatching", "boolval", MetaCommand::FLAG, false );

  command.SetOption ( "AlgorithmVerbosity", "d", false, "Algorithm verbosity (debug level)" );
  command.SetOptionLongTag ( "AlgorithmVerbosity", "verbose" );
  command.AddOptionField ( "AlgorithmVerbosity", "intval", MetaCommand::INT, false, "1" );
  command.SetOptionRange ( "AlgorithmVerbosity", "intval", "0", "100" );


  // Actually parse the command line

  if ( !command.Parse ( argc, argv ) )
  {
    exit ( EXIT_FAILURE );
  }


  // Store the parsed information into a struct
  args.fixedImageFile = command.GetValueAsString ( "FixedImageFile", "filename" );

  args.movingImageFile = command.GetValueAsString ( "MovingImageFile", "filename" );

  args.fixedImageMaskFile = command.GetValueAsString ( "FixedImageMaskFile", "filename" );

  args.movingImageMaskFile = command.GetValueAsString ( "MovingImageMaskFile", "filename" );

  args.inputFieldFile = command.GetValueAsString ( "InputFieldFile", "filename" );

  args.inputTransformFile = command.GetValueAsString ( "InputTransformFile", "filename" );

  args.outputImageFile = command.GetValueAsString ( "OutputImageFile", "filename" );

  args.outputDeformationFieldFile = command.GetValueAsString ( "OutputDeformationFieldFile", "filename" );

  args.outputInverseDeformationFieldFile = command.GetValueAsString ( "OutputInverseDeformationFieldFile", "filename" );

  args.outputVelocityFieldFile = command.GetValueAsString ( "OutputVelocityFieldFile", "filename" );
  
  args.outputJacobianFile = command.GetValueAsString ( "OutputJacobianFile", "filename" );

  unsigned int pos = args.outputImageFile.rfind ( "." );

  args.trueFieldFile = command.GetValueAsString ( "TrueFieldFile", "filename" );

  args.numIterations = parseUIntVector (
                         command.GetValueAsString ( "NumberOfIterationsPerLevels", "uintvect" ) );

  if ( args.numIterations.empty() || args.numIterations.size() > 10 )
  {
    std::cout << "NumberOfIterationsPerLevels.uintvect.size() : Value ("
              << args.numIterations.size() << ") is not in the range [1,10]" << std::endl;
    exit ( EXIT_FAILURE );
  }

  args.sigmaVel = command.GetValueAsFloat ( "VelocityFieldSigma", "floatval" );

  args.sigmaUp = command.GetValueAsFloat ( "UpdateFieldSigma", "floatval" );
  args.maxStepLength = command.GetValueAsFloat ( "MaximumUpdateStepLength", "floatval" );
  args.updateRule = command.GetValueAsInt ( "UpdateRule", "type" );
  args.gradientType = command.GetValueAsInt ( "GradientType", "type" );
  args.NumberOfBCHApproximationTerms = command.GetValueAsInt ( "NumberOfBCHApproximationTerms", "intval" );
  args.useHistogramMatching = command.GetValueAsBool ( "UseHistogramMatching", "boolval" );

  args.verbosity = 0;

  if ( command.GetOptionWasSet ( "AlgorithmVerbosity" ) )
  {
    args.verbosity = command.GetValueAsInt ( "AlgorithmVerbosity", "intval" );
  }
}


class SimpleCommandProgressUpdate : public itk::Command
{
  public:
    typedef SimpleCommandProgressUpdate   Self;
    typedef itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer; 
    itkNewMacro( Self );
    
    std::vector<bool>  _progress;
  protected:
    SimpleCommandProgressUpdate():
			_progress(11,false)
		{};
		
  public:
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      const itk::ProcessObject * filter =  dynamic_cast< const itk::ProcessObject * >( object );
      if( ! itk::ProgressEvent().CheckEvent( &event ) )
      {
        return;
      }
      
      int prg=filter->GetProgress()*(_progress.size()-1);
	if(prg>=_progress.size()) prg=_progress.size()-1;
	
	if(!_progress[prg])
	{
		std::cout << prg*100.0/(_progress.size()-1) << "% " <<std::flush;
		_progress[prg]=true;
	} else {
		std::cout << "."<<std::flush;
	}
    }
};



//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
template < class TPixel = float, unsigned int VImageDimension = 3 >
 class CommandIterationUpdate : public itk::Command
{
  public:
    typedef  CommandIterationUpdate                         Self;
    typedef  itk::Command                                   Superclass;
    typedef  itk::SmartPointer<Self>                        Pointer;

    typedef itk::Image< TPixel, VImageDimension >           InternalImageType;
    typedef itk::Vector< TPixel, VImageDimension >          VectorPixelType;
    typedef itk::Image<  VectorPixelType, VImageDimension > VelocityFieldType;
    typedef itk::Image<  VectorPixelType, VImageDimension > DeformationFieldType;

    typedef itk::LogDomainDeformableRegistrationFilterM <
    InternalImageType,
    InternalImageType,
    VelocityFieldType,
    InternalImageType >                                     LogDomainDeformableRegistrationFilterMType;

    typedef itk::MultiResolutionLogDomainDeformableRegistrationM <
    InternalImageType, InternalImageType,
    VelocityFieldType, InternalImageType, TPixel >         MultiResRegistrationFilterType;

    typedef itk::DisplacementFieldJacobianDeterminantFilter <
    DeformationFieldType, TPixel > JacobianFilterType;

    typedef itk::MinimumMaximumImageCalculator <
    InternalImageType >                                   MinMaxFilterType;

    typedef itk::WarpHarmonicEnergyCalculator <
    DeformationFieldType >                                HarmonicEnergyCalculatorType;

    typedef itk::VectorCentralDifferenceImageFunction <
    DeformationFieldType >                                WarpGradientCalculatorType;

    typedef typename WarpGradientCalculatorType::OutputType WarpGradientType;

    itkNewMacro ( Self );


    void SetTrueField ( const DeformationFieldType * truefield )
    {
      m_TrueField = truefield;

      m_TrueWarpGradientCalculator = WarpGradientCalculatorType::New();
      m_TrueWarpGradientCalculator->SetInputImage ( m_TrueField );

      m_CompWarpGradientCalculator =  WarpGradientCalculatorType::New();
    }

    void Execute ( itk::Object *caller, const itk::EventObject & event )
    {
      Execute ( ( const itk::Object * ) caller, event );
    }

    void Execute ( const itk::Object * object, const itk::EventObject & event )
    {
      if ( ! ( itk::IterationEvent().CheckEvent ( &event ) ) )
      {
        return;
      }

      typename DeformationFieldType::ConstPointer deffield = 0;

      unsigned int iter = -1;
      double metricbefore = -1.0;

      if ( const LogDomainDeformableRegistrationFilterMType * filter =
             dynamic_cast< const LogDomainDeformableRegistrationFilterMType * > ( object ) )
      {
        iter = filter->GetElapsedIterations() - 1;
        metricbefore = filter->GetMetric();
        deffield = const_cast<LogDomainDeformableRegistrationFilterMType *>
                   ( filter )->GetDeformationField();
      }
      else if ( const MultiResRegistrationFilterType * multiresfilter =
                  dynamic_cast< const MultiResRegistrationFilterType * > ( object ) )
      {
        std::cout << "Finished Multi-resolution iteration :" << multiresfilter->GetCurrentLevel() - 1 << std::endl;
        std::cout << "==============================" << std::endl << std::endl;
      }
      else
      {
        return;
      }

      if ( deffield )
      {
        std::cout << iter << ": MSE " << metricbefore << " - ";

        double fieldDist = -1.0;
        double fieldGradDist = -1.0;
        double tmp;

        if ( m_TrueField )
        {
          typedef itk::ImageRegionConstIteratorWithIndex<DeformationFieldType>
          FieldIteratorType;
          FieldIteratorType currIter (
            deffield, deffield->GetLargestPossibleRegion() );
          FieldIteratorType trueIter (
            m_TrueField, deffield->GetLargestPossibleRegion() );

          m_CompWarpGradientCalculator->SetInputImage ( deffield );

          fieldDist = 0.0;
          fieldGradDist = 0.0;

          for ( currIter.GoToBegin(), trueIter.GoToBegin();
                ! currIter.IsAtEnd(); ++currIter, ++trueIter )
          {
            fieldDist += ( currIter.Value() - trueIter.Value() ).GetSquaredNorm();

            // No need to add Id matrix here as we do a substraction
            tmp = (
                    ( m_CompWarpGradientCalculator->EvaluateAtIndex ( currIter.GetIndex() )
                      - m_TrueWarpGradientCalculator->EvaluateAtIndex ( trueIter.GetIndex() )
                    ).GetVnlMatrix() ).frobenius_norm();
            fieldGradDist += tmp * tmp;
          }

          fieldDist = sqrt ( fieldDist / ( double ) (

                               deffield->GetLargestPossibleRegion().GetNumberOfPixels() ) );
          fieldGradDist = sqrt ( fieldGradDist / ( double ) (
                                   deffield->GetLargestPossibleRegion().GetNumberOfPixels() ) );

          std::cout << "d(.,true) " << fieldDist << " - ";
          std::cout << "d(.,Jac(true)) " << fieldGradDist << " - ";
        }

        m_HarmonicEnergyCalculator->SetImage ( deffield );

        m_HarmonicEnergyCalculator->Compute();
        const double harmonicEnergy= m_HarmonicEnergyCalculator->GetHarmonicEnergy();
        std::cout << "harmo. " << harmonicEnergy << " - ";


        m_JacobianFilter->SetInput ( deffield );
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

        const double jacBelowZeroPrc = static_cast<double> ( jacBelowZero )
                                       / static_cast<double> ( numPix );


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


        std::cout << "max|Jac| " << maxJac << " - "
                  << "min|Jac| " << minJac << " - "
                  << "ratio(|Jac|<=0) " << jacBelowZeroPrc << std::endl;


        if ( this->m_Fid.is_open() )
        {
          if ( ! m_headerwritten )
          {
            this->m_Fid << "Iteration"
            << ", MSE before"
            << ", Harmonic energy"
            << ", min|Jac|"
            << ", 0.2% |Jac|"
            << ", 01% |Jac|"
            << ", 99% |Jac|"
            << ", 99.8% |Jac|"
            << ", max|Jac|"
            << ", ratio(|Jac|<=0)";

            if ( m_TrueField )
            {
              this->m_Fid << ", dist(warp,true warp)"
              << ", dist(Jac,true Jac)";
            }

            this->m_Fid << std::endl;

            m_headerwritten = true;
          }

          this->m_Fid << iter

          << ", " << metricbefore
          << ", " << harmonicEnergy
          << ", " << minJac
          << ", " << Q002
          << ", " << Q01
          << ", " << Q99
          << ", " << Q998
          << ", " << maxJac
          << ", " << jacBelowZeroPrc;

          if ( m_TrueField )
          {
            this->m_Fid << ", " << fieldDist
            << ", " << fieldGradDist;
          }

          this->m_Fid << std::endl;
        }
      }
    }

  protected:
    CommandIterationUpdate() :
        m_Fid ( "metricvalues.csv" ),
        m_headerwritten ( false )
    {
      m_JacobianFilter = JacobianFilterType::New();
      m_JacobianFilter->SetUseImageSpacing ( true );
      m_JacobianFilter->ReleaseDataFlagOn();

      m_Minmaxfilter = MinMaxFilterType::New();

      m_HarmonicEnergyCalculator = HarmonicEnergyCalculatorType::New();

      m_TrueField = 0;
      m_TrueWarpGradientCalculator = 0;
      m_CompWarpGradientCalculator = 0;
    };

    ~CommandIterationUpdate()
    {
      this->m_Fid.close();
    }

  private:
    std::ofstream m_Fid;
    bool m_headerwritten;
    typename JacobianFilterType::Pointer m_JacobianFilter;
    typename MinMaxFilterType::Pointer m_Minmaxfilter;
    typename HarmonicEnergyCalculatorType::Pointer m_HarmonicEnergyCalculator;
    typename DeformationFieldType::ConstPointer m_TrueField;
    typename WarpGradientCalculatorType::Pointer m_TrueWarpGradientCalculator;
    typename WarpGradientCalculatorType::Pointer m_CompWarpGradientCalculator;
};


template <unsigned int Dimension>
void LogDomainDemonsRegistrationFunction ( arguments args )
{
  // Declare the types of the images (float or double only)
  typedef float                               PixelType;
  typedef itk::Image< PixelType, Dimension >  ImageType;

  typedef itk::Vector< PixelType, Dimension > VectorPixelType;
  typedef typename itk::Image < VectorPixelType, Dimension > VelocityFieldType;
  typedef typename itk::Image < VectorPixelType, Dimension > DeformationFieldType;


  // Images we use
  typename ImageType::Pointer fixedImage = 0;
  typename ImageType::Pointer movingImage = 0;
  typename ImageType::Pointer fixedMask = 0;
  typename ImageType::Pointer movingMask = 0;
  typename VelocityFieldType::Pointer inputVelField = 0;


  // Set up the file readers
  typedef itk::ImageFileReader< ImageType >         FixedImageReaderType;
  typedef itk::ImageFileReader< ImageType >         MovingImageReaderType;
  typedef itk::ImageFileReader< ImageType >         MaskReaderType;
  typedef itk::ImageFileReader< VelocityFieldType > VelocityFieldReaderType;
  typedef itk::TransformFileReader                  TransformReaderType;

  {//for mem allocations

    typename FixedImageReaderType::Pointer fixedImageReader
    = FixedImageReaderType::New();
    typename MovingImageReaderType::Pointer movingImageReader
    = MovingImageReaderType::New();

    typename MaskReaderType::Pointer fixedMaskImageReader
    = MaskReaderType::New();

    typename MaskReaderType::Pointer movingMaskImageReader
    = MaskReaderType::New();

    fixedImageReader->SetFileName ( args.fixedImageFile.c_str() );
    movingImageReader->SetFileName ( args.movingImageFile.c_str() );

    if ( !args.fixedImageMaskFile.empty() )
      fixedMaskImageReader->SetFileName ( args.fixedImageMaskFile.c_str() );

    if ( !args.movingImageMaskFile.empty() )
      movingMaskImageReader->SetFileName ( args.movingImageMaskFile.c_str() );

    // Update the reader
    try
    {
      fixedImageReader->Update();
      movingImageReader->Update();

      if ( !args.fixedImageMaskFile.empty() )
      {
        fixedMaskImageReader->Update();
        fixedMask = fixedMaskImageReader->GetOutput();
        fixedMask->DisconnectPipeline();
      }

      if ( !args.movingImageMaskFile.empty() )
      {
        movingMaskImageReader->Update();
        movingMask = movingMaskImageReader->GetOutput();
        movingMask->DisconnectPipeline();
      }
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Could not read one of the input images." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }

    if ( ! args.inputFieldFile.empty() )
    {
      // Set up the file readers
      typename VelocityFieldReaderType::Pointer fieldReader = VelocityFieldReaderType::New();
      fieldReader->SetFileName ( args.inputFieldFile.c_str() );

      // Update the reader

      try
      {
        fieldReader->Update();
      }
      catch ( itk::ExceptionObject& err )
      {
        std::cout << "Could not read the input field." << std::endl;
        std::cout << err << std::endl;
        exit ( EXIT_FAILURE );
      }

      inputVelField = fieldReader->GetOutput();

      inputVelField->DisconnectPipeline();
    }
    
    if ( ! args.inputTransformFile.empty() )
    {
      std::string dumb1,dumb2;
      if(parse_xfm_file_name(args.inputTransformFile,dumb1,dumb2) ) //we are deling with .XFM filename
      {
#ifdef HAVE_MINC1 
        if( args.verbosity>0 )
          std::cout<<"Reading Deformations from XFM file..."<<std::endl;
        
        typedef typename minc::XfmTransform<double,Dimension,Dimension>  XfmTransformType;

        typename XfmTransformType::Pointer transform=XfmTransformType::New();
        
        transform->OpenXfm(args.inputTransformFile.c_str());
        
        SimpleCommandProgressUpdate::Pointer observer=SimpleCommandProgressUpdate::New();
	
  
      #if ( ITK_VERSION_MAJOR > 3 ) 
        typedef typename itk::TransformToDisplacementFieldSource<DeformationFieldType,double> TransformToDeformationSource;
      #else  
        typedef typename itk::TransformToDeformationFieldSource<DeformationFieldType,double> TransformToDeformationSource;
      #endif
        
        typename TransformToDeformationSource::Pointer transform_to_def=TransformToDeformationSource::New();
        transform_to_def->SetOutputParametersFromImage ( fixedImageReader->GetOutput() );
        transform_to_def->SetTransform(transform);
        transform_to_def->SetNumberOfThreads(1);
        
        if(args.verbosity>0)
        {
          transform_to_def->AddObserver ( itk::IterationEvent(), observer );
        }
        transform_to_def->Update();
        
        typedef typename itk::DisplacementToVelocityFieldLogFilter<DeformationFieldType,VelocityFieldType> DisplacementToVelocityFilter;
        typename DisplacementToVelocityFilter::Pointer displacement_to_velocity=DisplacementToVelocityFilter::New();
        displacement_to_velocity->SetInput(transform_to_def->GetOutput());
        
        if(args.verbosity>0)
        {
          displacement_to_velocity->AddObserver ( itk::IterationEvent(), observer );
        }
        
        displacement_to_velocity->Update();
        
        inputVelField = displacement_to_velocity->GetOutput();
        inputVelField->DisconnectPipeline();
        
        if( args.verbosity>0 )
          std::cout<<"Done..."<<std::endl;
#else
          std::cerr<<"Warning, XFM reader is not implemented in ITK4 yet!"<<std::endl;
#endif
          
      } else {
        typedef typename TransformReaderType::TransformType BaseTransformType;
        typename TransformReaderType::Pointer transformReader= TransformReaderType::New();
        BaseTransformType* baseTrsf ( 0 );
        // Set up the transform reader
        transformReader->SetFileName ( args.inputTransformFile.c_str() );

        // Update the reader

        try
        {
          transformReader->Update();
        }
        catch ( itk::ExceptionObject& err )
        {
          std::cout << "Could not read the input transform." << std::endl;
          std::cout << err << std::endl;
          exit ( EXIT_FAILURE );
        }
        

        const typename TransformReaderType::TransformListType* trsflistptr
        = transformReader->GetTransformList();

        if ( trsflistptr->empty() )
        {
          std::cout << "Could not read the input transform." << std::endl;
          exit ( EXIT_FAILURE );
        }
        else if ( trsflistptr->size() > 1 )
        {
          std::cout << "The input transform file contains more than one transform." << std::endl;
          exit ( EXIT_FAILURE );
        }

        baseTrsf = trsflistptr->front();

        if ( !baseTrsf )
        {
          std::cout << "Could not read the input transform." << std::endl;
          exit ( EXIT_FAILURE );
        }
        // Set up the TransformToDeformationFieldFilter
        typedef itk::TransformToVelocityFieldSource  <VelocityFieldType>  FieldGeneratorType;

        typedef typename FieldGeneratorType::TransformType TransformType;

        TransformType* trsf = dynamic_cast<TransformType*> ( baseTrsf );

        if ( !trsf )
        {
          std::cout << "Could not cast input transform to a usable transform." << std::endl;
          exit ( EXIT_FAILURE );
        }

        typename FieldGeneratorType::Pointer fieldGenerator = FieldGeneratorType::New();

        fieldGenerator->SetTransform ( trsf );
        fieldGenerator->SetOutputParametersFromImage ( fixedImageReader->GetOutput() );

        // Update the fieldGenerator

        try
        {
          fieldGenerator->Update();
        }
        catch ( itk::ExceptionObject& err )
        {
          std::cout << "Could not generate the input field." << std::endl;
          std::cout << err << std::endl;
          exit ( EXIT_FAILURE );
        }

        inputVelField = fieldGenerator->GetOutput();

        inputVelField->DisconnectPipeline();
      }
    }


    if ( !args.useHistogramMatching )
    {
      fixedImage = fixedImageReader->GetOutput();
      fixedImage->DisconnectPipeline();
      movingImage = movingImageReader->GetOutput();
      movingImage->DisconnectPipeline();
    }
    else
    {
      // match intensities
      ///\todo use inputDefField if any to get a better guess?
      typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> MatchingFilterType;
			
      typename MatchingFilterType::Pointer matcher = MatchingFilterType::New();

      matcher->SetInput ( movingImageReader->GetOutput() );
      matcher->SetReferenceImage ( fixedImageReader->GetOutput() );

      matcher->SetNumberOfHistogramLevels ( 1024 );
      matcher->SetNumberOfMatchPoints ( 7 );
      matcher->ThresholdAtMeanIntensityOn();

      // Update the matcher

      try
      {
        matcher->Update();
      }
      catch ( itk::ExceptionObject& err )
      {
        std::cout << "Could not match the input images." << std::endl;
        std::cout << err << std::endl;
        exit ( EXIT_FAILURE );
      }

      movingImage = matcher->GetOutput();

      movingImage->DisconnectPipeline();

      fixedImage = fixedImageReader->GetOutput();
      fixedImage->DisconnectPipeline();
    }

  }//end for mem allocations


  // Set up the demons filter output deformation field
  typename DeformationFieldType::Pointer defField = 0;
  typename DeformationFieldType::Pointer invDefField = 0;
  typename VelocityFieldType::Pointer velField = 0;

  {//for mem allocations

    // Set up the demons filter
    typedef typename itk::LogDomainDeformableRegistrationFilterM
    < ImageType, ImageType, VelocityFieldType>   BaseRegistrationFilterType;
    typename BaseRegistrationFilterType::Pointer filter;

    switch ( args.updateRule )
    {
      case 0:
      {
        // exp(v) <- exp(v) o exp(u) (log-domain demons)
        typedef typename itk::LogDomainDemonsRegistrationFilterM
        < ImageType, ImageType, VelocityFieldType>
        ActualRegistrationFilterType;
        typedef typename ActualRegistrationFilterType::GradientType GradientType;

        typename ActualRegistrationFilterType::Pointer actualfilter
        = ActualRegistrationFilterType::New();

        actualfilter->SetMaximumUpdateStepLength ( args.maxStepLength );
        actualfilter->SetUseGradientType (
          static_cast<GradientType> ( args.gradientType ) );
        actualfilter->SetNumberOfBCHApproximationTerms ( args.NumberOfBCHApproximationTerms );
        filter = actualfilter;

        break;
      }

      case 1:
      {
        // exp(v) <- Symmetrized( exp(v) o exp(u) ) (symmetriclog-domain demons)
        typedef typename itk::SymmetricLogDomainDemonsRegistrationFilterM
        < ImageType, ImageType, VelocityFieldType>
        ActualRegistrationFilterType;
        typedef typename ActualRegistrationFilterType::GradientType GradientType;

        typename ActualRegistrationFilterType::Pointer actualfilter
        = ActualRegistrationFilterType::New();

        actualfilter->SetMaximumUpdateStepLength ( args.maxStepLength );
        actualfilter->SetUseGradientType (
          static_cast<GradientType> ( args.gradientType ) );
        actualfilter->SetNumberOfBCHApproximationTerms ( args.NumberOfBCHApproximationTerms );
        filter = actualfilter;

        break;
      }

      default:
      {
        std::cout << "Unsupported update rule." << std::endl;
        exit ( EXIT_FAILURE );
      }
    }

    if ( args.sigmaVel > 0.1 )
    {
      filter->SmoothVelocityFieldOn();
      filter->SetStandardDeviations ( args.sigmaVel );
    }
    else
    {
      filter->SmoothVelocityFieldOff();
    }

    if ( args.sigmaUp > 0.1 )
    {
      filter->SmoothUpdateFieldOn();
      filter->SetUpdateFieldStandardDeviations ( args.sigmaUp );
    }
    else
    {
      filter->SmoothUpdateFieldOff();
    }

    //filter->SetIntensityDifferenceThreshold( 0.001 );

    if ( args.verbosity > 0 )
    {
      // Create the Command observer and register it with the registration filter.
      typename CommandIterationUpdate<PixelType, Dimension>::Pointer observer =
        CommandIterationUpdate<PixelType, Dimension>::New();

      if ( ! args.trueFieldFile.empty() )
      {
        if ( args.numIterations.size() > 1 )
        {
          std::cout << "You cannot compare the results with a true field in a multiresolution setting yet." << std::endl;
          exit ( EXIT_FAILURE );
        }

        // Set up the file readers
        typedef itk::ImageFileReader< DeformationFieldType > DeformationFieldReaderType;

        typename DeformationFieldReaderType::Pointer fieldReader = DeformationFieldReaderType::New();

        fieldReader->SetFileName ( args.trueFieldFile.c_str() );

        // Update the reader
        try
        {
          fieldReader->Update();
        }
        catch ( itk::ExceptionObject& err )
        {
          std::cout << "Could not read the true field." << std::endl;
          std::cout << err << std::endl;
          exit ( EXIT_FAILURE );
        }

        observer->SetTrueField ( fieldReader->GetOutput() );
      }

      filter->AddObserver ( itk::IterationEvent(), observer );
    }

    // Set up the multi-resolution filter
    //TODO: replace image mask type
    typedef typename itk::MultiResolutionLogDomainDeformableRegistrationM <
    ImageType, ImageType, VelocityFieldType, ImageType, PixelType >   MultiResRegistrationFilterType;

    typename MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();

    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction <
    VelocityFieldType, double > FieldInterpolatorType;

    typename FieldInterpolatorType::Pointer VectorInterpolator =
      FieldInterpolatorType::New();

    multires->GetFieldExpander()->SetInterpolator ( VectorInterpolator );

    multires->SetRegistrationFilter ( filter );

    multires->SetNumberOfLevels ( args.numIterations.size() );

    multires->SetNumberOfIterations ( &args.numIterations[0] );

    multires->SetFixedImage ( fixedImage );

    multires->SetMovingImage ( movingImage );

    if ( fixedMask )
      multires->SetFixedImageMask ( fixedMask );

    if ( movingMask )
      multires->SetMovingImageMask ( movingMask );

    
    multires->SetArbitraryInitialVelocityField ( inputVelField );


    if ( args.verbosity > 0 )
    {
      // Create the Command observer and register it with the registration filter.
      typename CommandIterationUpdate<PixelType, Dimension>::Pointer multiresobserver =
        CommandIterationUpdate<PixelType, Dimension>::New();
      multires->AddObserver ( itk::IterationEvent(), multiresobserver );
    }

    // Compute the deformation field
    try
    {
      multires->UpdateLargestPossibleRegion();
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }


    // Get various outputs

    // Final deformation field
    defField = multires->GetDeformationField();

    defField->DisconnectPipeline();

    // Inverse final deformation field
    invDefField = multires->GetInverseDeformationField();

    invDefField->DisconnectPipeline();

    // Final velocity field
    velField =  multires->GetVelocityField();

    velField->DisconnectPipeline();

  }//end for mem allocations


  typedef itk::WarpImageFilter < ImageType, ImageType, DeformationFieldType >  WarperType;
  typedef PixelType                                OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter   < ImageType, OutputImageType >                  CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


  if ( !args.outputImageFile.empty() )
  {
    typename WarperType::Pointer      warper = WarperType::New();
    typename WriterType::Pointer      writer =  WriterType::New();
    typename CastFilterType::Pointer  caster =  CastFilterType::New();

    warper->SetInput ( movingImage );
    warper->SetOutputSpacing ( fixedImage->GetSpacing() );
    warper->SetOutputOrigin ( fixedImage->GetOrigin() );
    warper->SetOutputDirection ( fixedImage->GetDirection() );
    
#if ( ITK_VERSION_MAJOR > 3 ) 
    warper->SetDisplacementField ( defField );
#else    
    warper->SetDeformationField ( defField );
#endif

    // warp the result
    // Write warped image out to file
    writer->SetFileName ( args.outputImageFile.c_str() );
    caster->SetInput ( warper->GetOutput() );
#ifdef HAVE_MINC4ITK
    mincify ( caster->GetOutput(), minc_history, NC_SHORT );
#else
    mincify ( caster->GetOutput(), minc_history, typeid(short).name() );
#endif
    writer->SetInput ( caster->GetOutput() );
    writer->SetUseCompression ( true );

    try
    {
      writer->Update();
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }
  }


  // Write output deformation field
  if ( !args.outputDeformationFieldFile.empty() )
  {
    int pos = args.outputDeformationFieldFile.rfind ( "." );
    std::string output_def_field = args.outputDeformationFieldFile;
    std::string out_grid_base_name;
    
    if ( parse_xfm_file_name(args.outputDeformationFieldFile,output_def_field,out_grid_base_name ))
    {

        std::ofstream out_xfm ( args.outputDeformationFieldFile.c_str() );

        out_xfm << "MNI Transform File" << std::endl << std::endl
        << "Transform_Type = Linear;" << std::endl
        << "Linear_Transform =" << std::endl
        << "    1 0 0 0" << std::endl
        << "    0 1 0 0" << std::endl
        << "    0 0 1 0;" << std::endl
        << "Transform_Type = Grid_Transform;" << std::endl
        << "Displacement_Volume = " << out_grid_base_name.c_str() << ";" << std::endl;
    }

    // Write the deformation field as an image of vectors.
    // Note that the file format used for writing the deformation field must be
    // capable of representing multiple components per pixel. This is the case
    // for the MetaImage and VTK file formats for example.
    typedef itk::ImageFileWriter< DeformationFieldType > FieldWriterType;

    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

    fieldWriter->SetFileName ( output_def_field.c_str() );

#ifdef HAVE_MINC4ITK
    mincify ( defField, minc_history, NC_SHORT );
#else
    mincify ( defField, minc_history, typeid(short).name() );
#endif
    fieldWriter->SetInput ( defField );

    fieldWriter->SetUseCompression ( true );

    try
    {
      fieldWriter->Update();
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }
  }

  // Write output inverse deformation field
  if ( !args.outputInverseDeformationFieldFile.empty() )
  {
    // Write the inverse deformation field as an image of vectors.
    // Note that the file format used for writing the inverse deformation field must be
    // capable of representing multiple components per pixel. This is the case
    // for the MetaImage and VTK file formats for example.
    typedef itk::ImageFileWriter< DeformationFieldType > FieldWriterType;
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName ( args.outputInverseDeformationFieldFile.c_str() );
#ifdef HAVE_MINC4ITK
    mincify ( invDefField, minc_history, NC_SHORT );
#else    
    mincify ( invDefField, minc_history, typeid(short).name() );
#endif
    fieldWriter->SetInput ( invDefField );
    fieldWriter->SetUseCompression ( true );

    try
    {
      fieldWriter->Update();
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }
  }

  // Write output inverse deformation field
  if ( !args.outputVelocityFieldFile.empty() )
  {
    // Write the velocity field as an image of vectors.
    // Note that the file format used for writing the velocity field must be
    // capable of representing multiple components per pixel. This is the case
    // for the MetaImage and VTK file formats for example.
    typedef itk::ImageFileWriter< VelocityFieldType > FieldWriterType;
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName ( args.outputVelocityFieldFile.c_str() );
#ifdef HAVE_MINC4ITK
    mincify ( velField, minc_history, NC_SHORT );
#else    
    mincify ( velField, minc_history, typeid(short).name() );
#endif    
    fieldWriter->SetInput ( velField );
    fieldWriter->SetUseCompression ( true );

    try
    {
      fieldWriter->Update();
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }
  }


  // Create and write warped grid image
  if ( args.verbosity > 0 )
  {
    typedef itk::Image< unsigned char, Dimension > GridImageType;
    typename GridImageType::Pointer gridImage = GridImageType::New();
    gridImage->SetRegions ( movingImage->GetRequestedRegion() );
    gridImage->SetOrigin ( movingImage->GetOrigin() );
    gridImage->SetSpacing ( movingImage->GetSpacing() );
    gridImage->Allocate();
    gridImage->FillBuffer ( 0 );

    typedef itk::ImageRegionIteratorWithIndex<GridImageType> GridImageIteratorWithIndex;
    GridImageIteratorWithIndex itergrid = GridImageIteratorWithIndex (
                                            gridImage, gridImage->GetRequestedRegion() );

    const int gridspacing ( 8 );

    for ( itergrid.GoToBegin(); !itergrid.IsAtEnd(); ++itergrid )
    {
      itk::Index<Dimension> index = itergrid.GetIndex();

      if ( Dimension == 2 || Dimension == 3 )
      {
        // produce an xy grid for all z
        if ( ( index[0] % gridspacing ) == 0 ||
             ( index[1] % gridspacing ) == 0 )
        {
          itergrid.Set ( itk::NumericTraits<unsigned char>::max() );
        }
      }
      else
      {
        unsigned int numGridIntersect = 0;

        for ( unsigned int dim = 0; dim < Dimension; dim++ )
        {
          numGridIntersect += ( ( index[dim] % gridspacing ) == 0 );
        }

        if ( numGridIntersect >= ( Dimension - 1 ) )
        {
          itergrid.Set ( itk::NumericTraits<unsigned char>::max() );
        }
      }
    }

    typedef itk::WarpImageFilter < GridImageType, GridImageType, DeformationFieldType >  GridWarperType;
    typename GridWarperType::Pointer gridwarper = GridWarperType::New();
    gridwarper->SetInput ( gridImage );
    gridwarper->SetOutputSpacing ( fixedImage->GetSpacing() );
    gridwarper->SetOutputOrigin ( fixedImage->GetOrigin() );
    gridwarper->SetOutputDirection ( fixedImage->GetDirection() );
    
#if ( ITK_VERSION_MAJOR > 3 ) 
    gridwarper->SetDisplacementField ( defField );
#else    
    gridwarper->SetDeformationField ( defField );
#endif

    // Write warped grid to file
    typedef itk::ImageFileWriter< GridImageType >  GridWriterType;

    typename GridWriterType::Pointer      gridwriter =  GridWriterType::New();
    gridwriter->SetFileName ( "WarpedGridImage.mnc" );
#ifdef HAVE_MINC4ITK
    mincify ( gridwarper->GetOutput(), minc_history,NC_BYTE );
#else
    mincify ( gridwarper->GetOutput(), minc_history,typeid(unsigned char).name() );
#endif    
    gridwriter->SetInput ( gridwarper->GetOutput() );
    gridwriter->SetUseCompression ( true );

    try
    {
      gridwriter->Update();
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }
  }


  // Create and write forewardwarped grid image
  if ( args.verbosity > 0 )
  {
    typedef itk::Image< unsigned char, Dimension > GridImageType;
    typedef itk::GridForwardWarpImageFilter<DeformationFieldType, GridImageType> GridForwardWarperType;

    typename GridForwardWarperType::Pointer fwWarper = GridForwardWarperType::New();
    fwWarper->SetInput ( defField );
    fwWarper->SetForegroundValue ( itk::NumericTraits<unsigned char>::max() );

    // Write warped grid to file
    typedef itk::ImageFileWriter< GridImageType >  GridWriterType;

    typename GridWriterType::Pointer      gridwriter =  GridWriterType::New();
    gridwriter->SetFileName ( "ForwardWarpedGridImage.mnc" );
    fwWarper->Update();
#ifdef HAVE_MINC4ITK
    mincify ( fwWarper->GetOutput(), minc_history, NC_BYTE );
#else    
    mincify ( fwWarper->GetOutput(), minc_history, typeid(unsigned char).name() );
#endif 
    gridwriter->SetInput ( fwWarper->GetOutput() );
    gridwriter->SetUseCompression ( true );

    try
    {
      gridwriter->Update();
    }
    catch ( itk::ExceptionObject& err )
    {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit ( EXIT_FAILURE );
    }
  }


  // compute final metric
  if ( args.verbosity > 0 )
  {

    typename WarperType::Pointer warper = WarperType::New();

    warper->SetInput ( movingImage );
    warper->SetOutputSpacing ( fixedImage->GetSpacing() );
    warper->SetOutputOrigin ( fixedImage->GetOrigin() );
    warper->SetOutputDirection ( fixedImage->GetDirection() );
    
#if ( ITK_VERSION_MAJOR > 3 ) 
    warper->SetDisplacementField ( defField );
#else    
    warper->SetDeformationField ( defField );
#endif

    warper->Update();

    double finalSSD = 0.0;
    typedef itk::ImageRegionConstIterator<ImageType> ImageConstIterator;

    ImageConstIterator iterfix = ImageConstIterator (
                                   fixedImage,
                                   fixedImage->GetRequestedRegion() );

    ImageConstIterator itermovwarp = ImageConstIterator (
                                       warper->GetOutput(),
                                       fixedImage->GetRequestedRegion() );

    for ( iterfix.GoToBegin(), itermovwarp.GoToBegin(); !iterfix.IsAtEnd(); ++iterfix, ++itermovwarp )
    {
      finalSSD += vnl_math_sqr ( iterfix.Get() - itermovwarp.Get() );
    }

    const double finalMSE = finalSSD / static_cast<double> (
                              fixedImage->GetRequestedRegion().GetNumberOfPixels() );

    std::cout << "MSE fixed image vs. warped moving image: " << finalMSE << std::endl;
  }



  // Create and write jacobian of the deformation field
  if ( args.verbosity > 0 || !args.outputJacobianFile.empty() )
  {
    typedef itk::DisplacementFieldJacobianDeterminantFilter
    <DeformationFieldType, PixelType> JacobianFilterType;
    typename JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
    jacobianFilter->SetInput ( defField );
    jacobianFilter->SetUseImageSpacing ( true );

    typename WriterType::Pointer      writer =  WriterType::New();
    typename CastFilterType::Pointer  caster =  CastFilterType::New();

    if(!args.outputJacobianFile.empty())
    {
      writer->SetFileName ( args.outputJacobianFile.c_str() );
      caster->SetInput ( jacobianFilter->GetOutput() );
#ifdef HAVE_MINC4ITK
      mincify ( jacobianFilter->GetOutput(), minc_history, NC_SHORT );
#else
      mincify ( jacobianFilter->GetOutput(), minc_history, typeid(short).name() );
#endif      
      writer->SetInput ( caster->GetOutput() );
      writer->SetUseCompression ( true );

      try
      {
	writer->Update();
      }
      catch ( itk::ExceptionObject& err )
      {
	std::cout << "Unexpected error." << std::endl;
	std::cout << err << std::endl;
	exit ( EXIT_FAILURE );
      }
    }
    
    if(args.verbosity>0)
    {
      typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxFilterType;

      typename MinMaxFilterType::Pointer minmaxfilter = MinMaxFilterType::New();
      minmaxfilter->SetImage ( jacobianFilter->GetOutput() );
      minmaxfilter->Compute();
      std::cout << "Minimum of the determinant of the Jacobian of the warp: "
		<< minmaxfilter->GetMinimum() << std::endl;
      std::cout << "Maximum of the determinant of the Jacobian of the warp: "
		<< minmaxfilter->GetMaximum() << std::endl;
    }
  }

}


int main ( int argc, char *argv[] )
{
  minc_history = minc_timestamp(argc,argv);
  
  struct arguments args;
  parseOpts ( argc, argv, args );

  std::cout << "Starting demons registration with the following arguments:" << std::endl;
  std::cout << args << std::endl << std::endl;

#if ITK_VERSION_MAJOR < 4 
  //Add Minc support
  itk::RegisterMincIO();
#endif  


  // FIXME uncomment for debug only
  // itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);

  // Get the image dimension
  itk::ImageIOBase::Pointer imageIO;

  try
  {
    imageIO = itk::ImageIOFactory::CreateImageIO (
                args.fixedImageFile.c_str(), itk::ImageIOFactory::ReadMode );

    if ( imageIO )
    {
      imageIO->SetFileName ( args.fixedImageFile.c_str() );
      imageIO->ReadImageInformation();
    }
    else
    {
      std::cout << "Could not read the fixed image information." << std::endl;
      exit ( EXIT_FAILURE );
    }
  }
  catch ( itk::ExceptionObject& err )
  {
    std::cout << "Could not read the fixed image information." << std::endl;
    std::cout << err << std::endl;
    exit ( EXIT_FAILURE );
  }

  switch ( imageIO->GetNumberOfDimensions() )
  {
    case 2:
      LogDomainDemonsRegistrationFunction<2> ( args );
      break;
    case 3:
      LogDomainDemonsRegistrationFunction<3> ( args );
      break;
    default:
      std::cout << "Unsuported dimension" << std::endl;
      exit ( EXIT_FAILURE );
  }

  return EXIT_SUCCESS;
}
