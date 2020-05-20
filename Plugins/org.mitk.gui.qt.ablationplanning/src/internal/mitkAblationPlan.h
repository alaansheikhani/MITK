/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef MITKABLATIONPLAN_H_HEADER_INCLUDED_
#define MITKABLATIONPLAN_H_HEADER_INCLUDED_
#include <itkDataObject.h>
#include <mitkCommon.h>
#include <mitkImage.h>
#include "mitkAblationUtils.h"

namespace mitk
{
  /**Documentation
   * \brief AblationPlan
   */
  class AblationPlan : public itk::DataObject
  {
  public:
    mitkClassMacroItkParent(AblationPlan, itk::DataObject);
    itkFactorylessNewMacro(Self);
    itkSetMacro(SegmentationImage, mitk::Image::Pointer);
    itkGetConstMacro(SegmentationImage, mitk::Image::Pointer);
    itkSetMacro(ImageDimension, mitk::Vector3D);
    itkGetConstMacro(ImageDimension, mitk::Vector3D);
    itkSetMacro(ImageSpacing, mitk::Vector3D);
    itkGetConstMacro(ImageSpacing, mitk::Vector3D);
    bool AddAblationZone(AblationUtils::AblationZone newZone);
    int GetNumberOfZones();
    AblationUtils::AblationZone* GetAblationZone(int id);
    bool RemoveAblationZone(AblationUtils::AblationZone &zone);
    bool RemoveAblationZone(int id);
    void DetectAndRemoveNotNeededVolumes();
    /** @return Returns 1 if plan b is better this, 0 if it is equal and -1 if it is worse*/
    int CompareTo(AblationPlan::Pointer b);

  protected:
    AblationPlan();

    ~AblationPlan() override;

  private:
    mitk::Image::Pointer m_SegmentationImage;
    mitk::Vector3D m_ImageDimension;
    mitk::Vector3D m_ImageSpacing;
    std::vector<AblationUtils::AblationZone> m_AblationZones;
  };
} // namespace mitk
#endif /* MITKABLATIONPLAN_H_HEADER_INCLUDED_ */
