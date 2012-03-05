/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPDEDeformableRegistrationFunctionM.h,v $
  Language:  C++
  Date:      $Date: 2009-01-26 21:45:56 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPDEDeformableRegistrationFunctionM_h
#define __itkPDEDeformableRegistrationFunctionM_h

#include "itkFiniteDifferenceFunction.h"

namespace itk {

/** \class PDEDeformableRegistrationFunctionM
 *
 * This is an abstract base class for all PDE functions which drives a
 * deformable registration algorithm. It is used by
 * PDEDeformationRegistrationFilter subclasses to compute the
 * output deformation field which will map a moving image onto
 * a fixed image.
 *
 * This class is templated over the fixed image type, moving image type
 * and the deformation field type.
 *
 * \sa PDEDeformableRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template<class TFixedImage, class TMovingImage, class TDeformationField,class TMask=TFixedImage>
class ITK_EXPORT PDEDeformableRegistrationFunctionM :
            public FiniteDifferenceFunction<TDeformationField>
{
public:
    /** Standard class typedefs. */
    typedef PDEDeformableRegistrationFunctionM           Self;
    typedef FiniteDifferenceFunction<TDeformationField> Superclass;
    typedef SmartPointer<Self>                          Pointer;
    typedef SmartPointer<const Self>                    ConstPointer;

    /** Run-time type information (and related methods) */
    itkTypeMacro( PDEDeformableRegistrationFunctionM,
                  FiniteDifferenceFunction );

    /** MovingImage image type. */
    typedef TMovingImage                            MovingImageType;
    typedef typename MovingImageType::ConstPointer  MovingImagePointer;

    /** FixedImage image type. */
    typedef TFixedImage                             FixedImageType;
    typedef typename FixedImageType::ConstPointer   FixedImagePointer;

    /** FixedImage image type. */
    typedef TMask                                   MaskImageType;
    typedef typename TMask::ConstPointer            MaskImagePointer;


    /** Deformation field type. */
    typedef TDeformationField    DeformationFieldType;
    typedef typename DeformationFieldType::Pointer
    DeformationFieldTypePointer;

    /** Set the moving image.  */
    void SetMovingImage( const MovingImageType * ptr )
    {
        m_MovingImage = ptr;
    }

    /** Get the moving image. */
    const MovingImageType * GetMovingImage(void) const
    {
        return m_MovingImage;
    }

    /** Set the fixed image. */
    void SetFixedImage( const FixedImageType * ptr )
    {
        m_FixedImage = ptr;
    }

    /** Get the fixed image. */
    const FixedImageType * GetFixedImage(void) const
    {
        return m_FixedImage;
    }

    /** Set the moving image.  */
    void SetMovingImageMask( const MaskImageType * ptr )
    {
        m_MovingImageMask = ptr;
    }

    /** Get the moving image. */
    const MaskImageType * GetMovingImageMask(void) const
    {
        return m_MovingImage;
    }

    /** Set the fixed image. */
    void SetFixedImageMask( const MaskImageType * ptr )
    {
        m_FixedImageMask = ptr;
    }

    /** Get the fixed image. */
    const MaskImageType * GetFixedImageMask(void) const
    {
        return m_FixedImageMask;
    }


    /** Set the deformation field image. */
    void SetDeformationField(  DeformationFieldTypePointer ptr )
    {
        m_DeformationField = ptr;
    }

    /** Get the deformation field. This function should have been
     *  declared const. It is not for backward compatibility reasons. */
    DeformationFieldType * GetDeformationField(void)
    {
        return m_DeformationField;
    }

    void SetEnergy( double e) {
        m_Energy=e;
    }
    double GetEnergy( ) const {
        return m_Energy;
    }
    void SetGradientStep( double e) {
        m_GradientStep = e;
    }
    double GetGradientStep( ) const {
        return m_GradientStep;
    }
    void SetNormalizeGradient( bool e) {
        m_NormalizeGradient=e;
    }
    bool GetNormalizeGradient( ) const {
        return m_NormalizeGradient;
    }

protected:
    PDEDeformableRegistrationFunctionM()
    {
        m_MovingImage = NULL;
        m_FixedImage = NULL;

        m_MovingImageMask = NULL;
        m_FixedImageMask = NULL;

        m_DeformationField = NULL;
        m_Energy = 0.0;
        m_NormalizeGradient = true;
        m_GradientStep = 1.0;
    }

    ~PDEDeformableRegistrationFunctionM() {}

    void PrintSelf(std::ostream& os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "MovingImage: ";
        os << m_MovingImage.GetPointer() << std::endl;
        os << indent << "FixedImage: ";
        os << m_FixedImage.GetPointer() << std::endl;

        os << indent << "MovingImageMask: ";
        os << m_MovingImageMask.GetPointer() << std::endl;
        os << indent << "FixedImageMask: ";
        os << m_FixedImageMask.GetPointer() << std::endl;

    }

    /** The moving image. */
    MovingImagePointer              m_MovingImage;

    /** The fixed image. */
    FixedImagePointer               m_FixedImage;

    /** The moving image. */
    MaskImagePointer                m_MovingImageMask;

    /** The fixed image. */
    MaskImagePointer                m_FixedImageMask;



    /** The deformation field. */
    DeformationFieldTypePointer     m_DeformationField;

    mutable double                  m_Energy;
    bool                            m_NormalizeGradient;
    mutable double                  m_GradientStep;
private:
    PDEDeformableRegistrationFunctionM(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

};


} // end namespace itk


#endif
