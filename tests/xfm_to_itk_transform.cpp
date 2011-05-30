#include <iostream>
#include "minc_helpers.h"

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <itkTransformFactory.h>
#include <itkMatrixOffsetTransformBase.h>


typedef itk::TransformFileWriter                TransformWriterType;

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
    
    itk::Matrix<double,3,3> rot;
    itk::Vector<double,3> tran;
    
    minc::read_linear_xfm(argv[1],rot, tran );
    AffineTransformType_d_3_3::Pointer transform=AffineTransformType_d_3_3::New();
    transform->SetMatrix(rot);
    transform->SetOffset(tran);
    
    itk::TransformFactory<AffineTransformType_d_3_3>::RegisterTransform();
    
    TransformWriterType::Pointer writer=TransformWriterType::New();
    writer->SetInput(transform);
    writer->SetFileName(argv[2]);
    writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 2;
  }
  return 0;
}
