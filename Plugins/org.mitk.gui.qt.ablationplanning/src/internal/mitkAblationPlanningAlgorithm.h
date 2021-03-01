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

#ifndef MITKABLATIONPLANNINGALGORITHM_H_HEADER_INCLUDED_
#define MITKABLATIONPLANNINGALGORITHM_H_HEADER_INCLUDED_
#include <itkObject.h>
#include <mitkCommon.h>
#include <mitkImage.h>
#include "mitkAblationUtils.h"
#include "mitkAblationPlan.h"
#include <omp.h>

namespace mitk
{
  /**Documentation
   * \brief AblationPlanningAlgorithm
   */
  class AblationPlanningAlgorithm : public itk::Object
  {
  public:
    mitkClassMacroItkParent(AblationPlanningAlgorithm, itk::Object);
    itkFactorylessNewMacro(Self);
    void ComputePlanning();
    itkGetMacro(AblationPlan,mitk::AblationPlan::Pointer);
    void SetAdjustableParameters(int iterations, double maxAblationRadius, double ablationRadius, double minAblationRadius, double toleranceNonAblatedTumorSafetyMarginVolume);
    void SetSegmentationData(mitk::Image::Pointer segmentationImage, mitk::Vector3D imageDimension, mitk::Vector3D imageSpacing);
    void SetStartingPoint(mitk::Point3D startingPositionInWorldCoordinates, itk::Index<3> startingPositionIndexCoordinates);
    void SetSafetyMargin(std::vector<itk::Index<3>> tumorTissueSafetyMarginIndices);

  protected:
    AblationPlanningAlgorithm();

    ~AblationPlanningAlgorithm() override;

  private:
    /* Algorithm adjustable Parameters */
    int m_Iterations;
    double m_MaxAblationRadius; // Maximal ablation radius
    double m_AblationRadius; // Desired ablation radius
    double m_MinAblationRadius; // Minimal ablation radius
    double m_ToleranceNonAblatedTumorSafetyMarginVolume; // Tolerance Non-Ablated Volume

    /* Segmentation data */
    mitk::Image::Pointer m_SegmentationImage;
    mitk::Vector3D m_ImageDimension;
    mitk::Vector3D m_ImageSpacing;

    /* Safety margin */
     /*!
     * \brief Vector storing the index coordinates of all pixels, which are tumor tissue or
     * safety margin.
     */
    std::vector<itk::Index<3>> m_TumorTissueSafetyMarginIndices;

    /* Starting point */
    //mitk::Point3D m_AblationStartingPositionInWorldCoordinates; //wozu?
    //itk::Index<3> m_AblationStartingPositionIndexCoordinates; //wozu?
    mitk::Point3D m_TempAblationStartingPositionInWorldCoordinates;
    itk::Index<3> m_TempAblationStartingPositionIndexCoordinates;
    bool m_ManualAblationStartingPositionSet;

    /* Output / result */
    mitk::AblationPlan::Pointer m_AblationPlan;


    /*!
     * \brief Temporary vector storing the index coordinates of all circle centers of the ablation zones
     * when calculating the best ablation zone distribution.
     */
    std::vector<mitk::AblationZone> m_TempAblationZones;
  };
} // namespace mitk
#endif /* MITKABLATIONPLANNINGALGORITHM_H_HEADER_INCLUDED_ */
