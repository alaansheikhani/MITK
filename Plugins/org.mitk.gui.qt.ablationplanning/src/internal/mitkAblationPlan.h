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
#include "mitkAblationZone.h"
#include <vector>

namespace mitk
{
  /**Documentation
   * \brief AblationPlan
   */
  class AblationPlan : public itk::DataObject
  {
  public:
    struct AblationPlanStatistics{
      int tumorVolume;
      int safetyMarginVolume;
      int tumorAndSafetyMarginVolume;
      int totalAblationVolume;
      int ablationVolumeAblatedMoreThanOneTime;
      double factorOverlappingAblationZones;
      double factorAblatedVolumeOutsideSafetyMargin;
      double factorNonAblatedVolume;
      double factorMaxNonAblatedVolume;
      std::vector<std::vector<itk::Index<3>>> agglomerations;
    };  
    AblationPlanStatistics GetStatistics();
    void SetStatistics(AblationPlanStatistics s);
    itkGetConstMacro(StatsSet, bool);
    mitkClassMacroItkParent(AblationPlan, itk::DataObject);
    itkFactorylessNewMacro(Self);
    void SetSegmentationImage(mitk::Image::Pointer s);
    mitk::Image::Pointer GetSegmentationImage();
    itkSetMacro(ImageDimension, mitk::Vector3D);
    itkGetConstMacro(ImageDimension, mitk::Vector3D);
    itkSetMacro(ImageSpacing, mitk::Vector3D);
    itkGetConstMacro(ImageSpacing, mitk::Vector3D);
    bool AddAblationZone(mitk::AblationZone newZone);
    int GetNumberOfZones();
    mitk::AblationZone* GetAblationZone(int id);
    bool RemoveAblationZone(mitk::AblationZone &zone);
    bool RemoveAblationZone(int id);
    /** @return Returns 1 if plan b is better this, 0 if it is equal and -1 if it is worse*/
    int CompareTo(AblationPlan::Pointer b);
    void CalculcateSolutionValue();
    double CalculateSolutionValueOfZoneNumbers();
    //@returns Returns between 1 (lowest number of Ablation Zones) and 0 (higest number of Ablation Zones)
    double CalculateSolutionValueOfOverlappingZones();
    double CalculateSolutionValueOfZoneDifference();
    double CalculateSolutionValueOfVolumeOutsideFactor();
    double CalculateSolutionValueOfNonAblated();
    void SetSolutionValue(double add);
    double GetSolutionValue();
    void SetMaxAblationZoneNumber(int n);
    void SetMinAblationZoneNumber(int n);
    int GetMinAblationZoneNumber();
    int GetMaxAblationZoneNumber();
    void SetMaxOverlappingZonesFactor(double overlapFactor);
    void SetMinOverlappingZonesFactor(double overlapFactor);
    double GetMaxOverlappingZonesFactor();
    double GetMinOverlappingZonesFactor();
    void SetMaxVolumeOutsideFactor(double factorOutside);
    void SetMinVolumeOutsideFactor(double factorOutside);
    double GetMaxVolumeOutsideFactor();
    double GetMinVolumeOutsideFactor();

  protected:
    AblationPlan();

    ~AblationPlan() override;

  private:
    mitk::Image::Pointer m_SegmentationImage;
    mitk::Vector3D m_ImageDimension;
    mitk::Vector3D m_ImageSpacing;
    std::vector<mitk::AblationZone> m_AblationZones;
    bool m_StatsSet;
    AblationPlanStatistics m_Stats;
    double m_SolutionValue;
    int m_MaxAblationZoneNumber;
    int m_MinAblationZoneNumber;
    double m_MaxOverlappingZonesFactor;
    double m_MinOverlappingZonesFactor;
    double m_MaxVolumeOutsideFactor;
    double m_MinVolumeOutsideFactor;
  };
} // namespace mitk
#endif /* MITKABLATIONPLAN_H_HEADER_INCLUDED_ */
