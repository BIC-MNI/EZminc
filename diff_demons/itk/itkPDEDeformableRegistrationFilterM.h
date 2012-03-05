/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPDEDeformableRegistrationFilterM.h,v $
  Language:  C++
  Date:      $Date: 2009-10-19 15:18:52 $
  Version:   $Revision: 1.26 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPDEDeformableRegistrationFilterM_h
#define __itkPDEDeformableRegistrationFilterM_h

#include "itkDenseFiniteDifferenceImageFilter.h"
#include "itkPDEDeformableRegistrationFunctionM.h"

namespace itk {

    /**
     * \class PDEDeformableRegistrationFilterM
     * \brief Deformably register two images using a PDE algorithm.
     *
     * PDEDeformableRegistrationFilterM is a base case for filter implementing
     * a PDE deformable algorithm that register two images by computing the
     * deformation field which will map a moving image onto a fixed image.
     *
     * A deformation field is represented as a image whose pixel type is some
     * vector type with at least N elements, where N is the dimension of
     * the fixed image. The vector type must support element access via operator
     * []. It is assumed that the vector elements behave like floating point
     * scalars.
     *
     * This class is templated over the fixed image type, moving image type
     * and the deformation Field type.
     *
     * The input fixed and moving images are set via methods SetFixedImage
     * and SetMovingImage respectively. An initial deformation field maybe set via
     * SetInitialDeformationField or SetInput. If no initial field is set,
     * a zero field is used as the initial condition.
     *
     * The output deformation field can be obtained via methods GetOutput
     * or GetDeformationField.
     *
     * The PDE algorithm is run for a user defined number of iterations.
     * Typically the PDE algorithm requires period Gaussin smoothing of the
     * deformation field to enforce an elastic-like condition. The amount
     * of smoothing is governed by a set of user defined standard deviations
     * (one for each dimension).
     *
     * In terms of memory, this filter keeps two internal buffers: one for storing
     * the intermediate updates to the field and one for double-buffering when
     * smoothing the deformation field. Both buffers are the same type and size as the
     * output deformation field.
     *
     * This class make use of the finite difference solver hierarchy. Update
     * for each iteration is computed using a PDEDeformableRegistrationFunction.
     *
     * \warning This filter assumes that the fixed image type, moving image type
     * and deformation field type all have the same number of dimensions.
     *
     * \sa PDEDeformableRegistrationFunction.
     * \ingroup DeformableImageRegistration
     */
    template<class TFixedImage, class TMovingImage, class TDeformationField,class TMask=TFixedImage>
    class ITK_EXPORT PDEDeformableRegistrationFilterM :
                public DenseFiniteDifferenceImageFilter<TDeformationField,TDeformationField> {
        public:
            /** Standard class typedefs. */
            typedef PDEDeformableRegistrationFilterM    Self;
            typedef DenseFiniteDifferenceImageFilter<
            TDeformationField,TDeformationField>     Superclass;
            typedef SmartPointer<Self>                 Pointer;
            typedef SmartPointer<const Self>           ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro ( Self );

            /** Run-time type information (and related methods) */
            itkTypeMacro ( PDEDeformableRegistrationFilterM,
                           DenseFiniteDifferenceImageFilter );

            /** FixedImage image type. */
            typedef TFixedImage                            FixedImageType;
            typedef typename FixedImageType::Pointer       FixedImagePointer;
            typedef typename FixedImageType::ConstPointer  FixedImageConstPointer;

            /** MovingImage image type. */
            typedef TMovingImage                           MovingImageType;
            typedef typename MovingImageType::Pointer      MovingImagePointer;
            typedef typename MovingImageType::ConstPointer MovingImageConstPointer;

            /** Deformation field type. */
            typedef TDeformationField                       DeformationFieldType;
            typedef typename DeformationFieldType::Pointer  DeformationFieldPointer;

            /** Mask image type. */
            typedef TMask                                   MaskImageType;
            typedef typename TMask::ConstPointer            MaskImagePointer;


            /** Types inherithed from the superclass */
            typedef typename Superclass::OutputImageType    OutputImageType;

            /** FiniteDifferenceFunction type. */
            typedef typename Superclass::FiniteDifferenceFunctionType
            FiniteDifferenceFunctionType;

            /** PDEDeformableRegistrationFilterMFunction type. */
            typedef PDEDeformableRegistrationFunctionM<FixedImageType,MovingImageType,
            DeformationFieldType,MaskImageType>  PDEDeformableRegistrationFunctionType;

            /** Inherit some enums and typedefs from the superclass. */
            itkStaticConstMacro ( ImageDimension, unsigned int,
                                  Superclass::ImageDimension );

            /** Set the fixed image. */
            void SetFixedImage ( const FixedImageType * ptr );

            /** Get the fixed image. */
            const FixedImageType * GetFixedImage ( void ) const;

            /** Set the moving image. */
            void SetMovingImage ( const MovingImageType * ptr );

            /** Get the moving image. */
            const MovingImageType * GetMovingImage ( void ) const;


            /** Set the moving image.  */
            void SetMovingImageMask ( const MaskImageType * ptr ) {
                m_MovingImageMask = ptr;
            }

            /** Get the moving image. */
            const MaskImageType * GetMovingImageMask ( void ) const {
                return m_MovingImageMask;
            }

            /** Set the fixed image. */
            void SetFixedImageMask ( const MaskImageType * ptr ) {
                m_FixedImageMask = ptr;
            }

            /** Get the fixed image. */
            const MaskImageType * GetFixedImageMask ( void ) const {
                return m_FixedImageMask;
            }


            /** Set initial deformation field. */
            void SetInitialDeformationField ( const DeformationFieldType * ptr ) {
                this->SetInput ( ptr );
            }

            /** Get output deformation field. */
            DeformationFieldType * GetDeformationField() {
                return this->GetOutput();
            }

            /** Get the number of valid inputs.  For PDEDeformableRegistration,
             * this checks whether the fixed and moving images have been
             * set. While PDEDeformableRegistration can take a third input as an
             * initial deformation field, this input is not a required input.
             */
            virtual std::vector<SmartPointer<DataObject> >::size_type GetNumberOfValidRequiredInputs() const;

            /** Set/Get whether the deformation field is smoothed
             * (regularized). Smoothing the deformation yields a solution
             * elastic in nature. If SmoothDeformationField is on, then the
             * deformation field is smoothed with a Gaussian whose standard
             * deviations are specified with SetStandardDeviations() */
            itkSetMacro ( SmoothDeformationField, bool );
            itkGetConstMacro ( SmoothDeformationField, bool );
            itkBooleanMacro ( SmoothDeformationField );

            /** Set the Gaussian smoothing standard deviations for the
             * deformation field. The values are set with respect to pixel
             * coordinates. */
            itkSetVectorMacro ( StandardDeviations, double, ImageDimension );
            virtual void SetStandardDeviations ( double value );

            /** Get the Gaussian smoothing standard deviations use for smoothing
             * the deformation field. */
            const double * GetStandardDeviations ( void ) {
                return ( double * ) m_StandardDeviations;
            }

            /** Set/Get whether the update field is smoothed
             * (regularized). Smoothing the update field yields a solution
             * viscous in nature. If SmoothUpdateField is on, then the
             * update field is smoothed with a Gaussian whose standard
             * deviations are specified with SetUpdateFieldStandardDeviations() */
            itkSetMacro ( SmoothUpdateField, bool );
            itkGetConstMacro ( SmoothUpdateField, bool );
            itkBooleanMacro ( SmoothUpdateField );

            /** Set the Gaussian smoothing standard deviations for the update
             * field. The values are set with respect to pixel coordinates. */
            itkSetVectorMacro ( UpdateFieldStandardDeviations, double, ImageDimension );
            virtual void SetUpdateFieldStandardDeviations ( double value );

            /** Get the Gaussian smoothing standard deviations used for
             * smoothing the update field. */
            const double * GetUpdateFieldStandardDeviations ( void ) {
                return ( double * ) m_UpdateFieldStandardDeviations;
            }



            /** Stop the registration after the current iteration. */
            virtual void StopRegistration() {
                m_StopRegistrationFlag = true;
            }

            /** Set/Get the desired maximum error of the Guassian kernel approximate.
             * \sa GaussianOperator. */
            itkSetMacro ( MaximumError, double );
            itkGetConstMacro ( MaximumError, double );

            /** Set/Get the desired limits of the Gaussian kernel width.
             * \sa GaussianOperator. */
            itkSetMacro ( MaximumKernelWidth, unsigned int );
            itkGetConstMacro ( MaximumKernelWidth, unsigned int );


        protected:
            PDEDeformableRegistrationFilterM();
            ~PDEDeformableRegistrationFilterM() {}
            void PrintSelf ( std::ostream& os, Indent indent ) const;

            /** Supplies the halting criteria for this class of filters.  The
             * algorithm will stop after a user-specified number of iterations. */
            virtual bool Halt() {

                if ( m_StopRegistrationFlag ) {
                    return true;
                }

                return this->Superclass::Halt();
            }

            /** A simple method to copy the data from the input to the output.
             * If the input does not exist, a zero field is written to the output. */
            virtual void CopyInputToOutput();

            /** Initialize the state of filter and equation before each iteration.
             * Progress feeback is implemented as part of this method. */
            virtual void InitializeIteration();

            /** Utility to smooth the deformation field (represented in the Output)
             * using a Guassian operator. The amount of smoothing can be specified
             * by setting the StandardDeviations. */
            //virtual void SmoothDisplacementField();
            virtual void SmoothDeformationField();

            /** Utility to smooth the UpdateBuffer using a Gaussian operator.
             * The amount of smoothing can be specified by setting the
             * UpdateFieldStandardDeviations. */
            virtual void SmoothUpdateField();

            /** This method is called after the solution has been generated. In this case,
             * the filter release the memory of the internal buffers. */
            virtual void PostProcessOutput();

            /** This method is called before iterating the solution. */
            virtual void Initialize();

            /** By default the output deformation field has the same Spacing, Origin
             * and LargestPossibleRegion as the input/initial deformation field.  If
             * the initial deformation field is not set, the output information is
             * copied from the fixed image. */
            virtual void GenerateOutputInformation();

            /** It is difficult to compute in advance the input moving image region
             * required to compute the requested output region. Thus the safest
             * thing to do is to request for the whole moving image.
             *
             * For the fixed image and deformation field, the input requested region
             * set to be the same as that of the output requested region. */
            virtual void GenerateInputRequestedRegion();

        private:
            PDEDeformableRegistrationFilterM ( const Self& ); //purposely not implemented
            void operator= ( const Self& ); //purposely not implemented

            /** Standard deviation for Gaussian smoothing */
            double m_StandardDeviations[ImageDimension];
            double m_UpdateFieldStandardDeviations[ImageDimension];

            /** Modes to control smoothing of the update and deformation fields */
            bool m_SmoothDeformationField;
            bool m_SmoothUpdateField;

            /** Temporary deformation field use for smoothing the
             * the deformation field. */
            DeformationFieldPointer   m_TempField;

        private:
            /** Maximum error for Gaussian operator approximation. */
            double                    m_MaximumError;

            /** Limits of Guassian kernel width. */
            unsigned int              m_MaximumKernelWidth;

            /** Flag to indicate user stop registration request. */
            bool                      m_StopRegistrationFlag;

            /** The moving image. */
            MaskImagePointer                m_MovingImageMask;

            /** The fixed image. */
            MaskImagePointer                m_FixedImageMask;

    };


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPDEDeformableRegistrationFilterM.txx"
#endif

#endif
