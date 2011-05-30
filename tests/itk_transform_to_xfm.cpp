#include <iostream>
#include "minc_helpers.h"
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <itkTransformFactory.h>
#include <itkMatrixOffsetTransformBase.h>


typedef itk::TransformFileReader                TransformReaderType;

typedef itk::MatrixOffsetTransformBase<double, 3,3> AffineTransformType_d_3_3;
typedef itk::MatrixOffsetTransformBase<float, 3,3>  AffineTransformType_f_3_3;

int main(int argc,char **argv)
{
  if(argc<2)
  {
    std::cerr<<"Usage:"<<argv[0]<<" <in.txt> <out.xfm>"<<std::endl;
    return 1;
  }
  
  try
  {
    itk::TransformFactory<AffineTransformType_f_3_3>::RegisterTransform();
    itk::TransformFactory<AffineTransformType_d_3_3>::RegisterTransform();
    
    TransformReaderType::Pointer reader=TransformReaderType::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    
    AffineTransformType_d_3_3 * affine_d;
    AffineTransformType_f_3_3 * affine_f;
    

    if( (affine_d = dynamic_cast< AffineTransformType_d_3_3 * >( reader->GetTransformList()->front().GetPointer() )) != NULL )
    {
      std::cout << "It is MatrixOffsetTransformBase<double, 3,3> " << std::endl;
      minc::write_linear_xfm(argv[2],affine_d->GetMatrix(),affine_d->GetTranslation());
      
    } else {
      std::cerr << "Transformation format usupported at the moment " << std::endl;
    }

    
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
  return 0;
}
