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

template< class TPixel = float, unsigned int VImageDimension = 3 >
class CommandIterationUpdate:public itk::Command
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

  itkNewMacro(Self);

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event );
  }

  void Execute(const itk::Object *object, const itk::EventObject & event)
  {
    if ( !( itk::IterationEvent().CheckEvent(&event) ) )
      {
      return;
      }

    if ( VelocitorFilterType * filter = dynamic_cast< VelocitorFilterType * >(
            const_cast< Object * >( object ) ) )
      {
      typename VelocityFieldType::Pointer velField = filter->GetOutput();

      typename ExponentiatorFilterType::Pointer exponentiator = ExponentiatorFilterType::New();
      exponentiator->SetInput (velField);

      try
        {
        exponentiator->Update();
        }
      catch ( itk::ExceptionObject & e )
        {
        std::cerr << e;
        itkExceptionMacro (<< "Error in exponentiation");
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

int main(int argc, char *argv[])
{
  if ( argc < 3 )
    {
    std::cerr << "Usage:" << std::endl;
    std::cerr << argv[0] << " input output number_of_iterations" << std::endl;
    return -1;
    }

  typedef itk::Vector< double, 3 >    VectorType;
  typedef itk::Image< VectorType, 3 > ImageType;

  typedef itk::ImageFileReader< ImageType >                                   ReaderType;
  typedef itk::ImageFileWriter< ImageType >                                   WriterType;
  typedef itk::DisplacementToVelocityFieldLogFilter< ImageType, ImageType >   VelocitorFilterType;
  typedef itk::ExponentialDeformationFieldImageFilter< ImageType, ImageType > ExponentiatorFilterType;
  typedef itk::DisplacementFieldCompositionFilter< ImageType, ImageType >     CompositionFilterType;

  std::cout << "Reading...";
  ImageType::Pointer displacement = 0;
    {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName (argv[1]);
    try
      {
      reader->Update();
      }
    catch ( itk::ExceptionObject & e )
      {
      std::cerr << e;
      return -1;
      }

    displacement = reader->GetOutput();
    displacement->DisconnectPipeline();
    }
  std::cout << "Done." << std::endl;

  std::cout << "Velocitying...";
  ImageType::Pointer velocity = 0;
    {
    VelocitorFilterType::Pointer velocitor = VelocitorFilterType::New();
    velocitor->SetInput (displacement);
    velocitor->SetNumberOfIterations ( std::stoi(argv[3]) );
    velocitor->SmoothVelocityFieldOn();
    velocitor->SetSigma (2.0);
    velocitor->SetNumberOfExponentialIntegrationSteps (500);
    velocitor->SetNumberOfBCHApproximationTerms (3);

    CommandIterationUpdate< double, 3 >::Pointer observer = CommandIterationUpdate< double, 3 >::New();
    velocitor->AddObserver (itk::IterationEvent(), observer);

    try
      {
      velocitor->Update();
      }
    catch ( itk::ExceptionObject & e )
      {
      std::cerr << e;
      return -1;
      }

    velocity = velocitor->GetOutput();
    velocity->DisconnectPipeline();
    }
  std::cout << "Done." << std::endl;

  std::cout << "Inversing...";
  ImageType::Pointer inverse = 0;
    {
    ExponentiatorFilterType::Pointer exponentiator = ExponentiatorFilterType::New();
    exponentiator->SetInput (velocity);
    exponentiator->ComputeInverseOn();

    try
      {
      exponentiator->Update();
      }
    catch ( itk::ExceptionObject & e )
      {
      std::cerr << e;
      return -1;
      }

    inverse = exponentiator->GetOutput();
    inverse->DisconnectPipeline();
    }
  std::cout << "Done." << std::endl;

  std::cout << "Saving...";
    {
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName (argv[2]);
    writer->SetInput (inverse);
    try
      {
      writer->Update();
      }
    catch ( itk::ExceptionObject & e )
      {
      std::cerr << e;
      return -1;
      }
    }
  std::cout << "Done." << std::endl;

  return EXIT_SUCCESS;
}
