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

#include "mitkAblationPlanningAlgorithm.h"

mitk::AblationPlanningAlgorithm::AblationPlanningAlgorithm(){}

mitk::AblationPlanningAlgorithm::~AblationPlanningAlgorithm() {}

void mitk::AblationPlanningAlgorithm::SetAdjustableParameters(int iterations,double maxAblationRadius, double ablationRadius, double minAblationRadius, double toleranceNonAblatedTumorSafetyMarginVolume){
  m_Iterations = iterations;
  m_MaxAblationRadius = maxAblationRadius;
  m_AblationRadius = ablationRadius;
  m_MinAblationRadius = minAblationRadius;
  m_ToleranceNonAblatedTumorSafetyMarginVolume = toleranceNonAblatedTumorSafetyMarginVolume;
}

void mitk::AblationPlanningAlgorithm::SetSegmentationData(mitk::Image::Pointer segmentationImage, mitk::Vector3D imageDimension, mitk::Vector3D imageSpacing){
  m_SegmentationImage = segmentationImage;
  m_ImageDimension = imageDimension;
  m_ImageSpacing = imageSpacing;
}

void mitk::AblationPlanningAlgorithm::SetStartingPoint(mitk::Point3D startingPositionInWorldCoordinates, itk::Index<3> startingPositionIndexCoordinates){
  m_TempAblationStartingPositionInWorldCoordinates = startingPositionInWorldCoordinates;
  m_TempAblationStartingPositionIndexCoordinates = startingPositionIndexCoordinates;
  m_ManualAblationStartingPositionSet = true;
}

void mitk::AblationPlanningAlgorithm::SetSafetyMargin(std::vector<itk::Index<3>> tumorTissueSafetyMarginIndices){
  m_TumorTissueSafetyMarginIndices = tumorTissueSafetyMarginIndices;
}

void mitk::AblationPlanningAlgorithm::ComputePlanning(){
//==================== Find ablation proposal by iteratively create and test random propsals
  //==========================================
  AblationUtils::FillVectorContainingIndicesOfTumorTissueSafetyMargin(
    m_SegmentationImage, m_ImageDimension, m_TumorTissueSafetyMarginIndices);

  std::vector<mitk::AblationPlan::Pointer> AllFoundPlans = std::vector<mitk::AblationPlan::Pointer>();
  //AllFoundPlans.resize(m_Controls.repititionsCalculatingAblationZonesSpinBox->value());

  // Start of for-loop (main iterative loop, each iteration = one proposal):
  for (int iteration = 1; iteration <= m_Iterations; ++iteration)
  {
    mitk::AblationPlan::Pointer currentPlan = mitk::AblationPlan::New();
    currentPlan->SetSegmentationImage(m_SegmentationImage->Clone());
    currentPlan->SetImageDimension(m_ImageDimension);
    currentPlan->SetImageSpacing(m_ImageSpacing);

    double startingZoneRadius = 0;
    MITK_INFO << "Iteration: " << iteration;
    if (!m_ManualAblationStartingPositionSet || iteration > 1)
    {
      QString position = AblationUtils::FindAblationStartingPosition(currentPlan->GetSegmentationImage(),
                                                                     m_TumorTissueSafetyMarginIndices,
                                                                     m_AblationRadius,
                                                                     m_MaxAblationRadius,
                                                                     m_TempAblationStartingPositionIndexCoordinates,
                                                                     m_TempAblationStartingPositionInWorldCoordinates,
                                                                     startingZoneRadius,
                                                                     m_ImageDimension,
                                                                     m_ImageSpacing);
      //m_Controls.ablationStartingPointLabel->setText(position);
      MITK_INFO << "Found starting point: " << position.toLatin1().toStdString();
    }
    else
    {
      startingZoneRadius = this->m_AblationRadius;
    }

    //------------ Random distribution model calculations: ---------------
    {
      double size = m_TumorTissueSafetyMarginIndices.size();
      std::vector<itk::Index<3>> indices = m_TumorTissueSafetyMarginIndices;
      AblationUtils::CalculateAblationVolume(m_TempAblationStartingPositionIndexCoordinates,
                                             currentPlan->GetSegmentationImage(),
                                             startingZoneRadius,
                                             m_ImageSpacing,
                                             m_ImageDimension,
                                             m_TempAblationZones);
      AblationUtils::RemoveAblatedPixelsFromGivenVector(m_TempAblationZones.at(0).indexCenter,
                                                        indices,
                                                        currentPlan->GetSegmentationImage(),
                                                        m_TempAblationZones.at(0).radius,
                                                        m_ImageDimension,
                                                        m_ImageSpacing);
      currentPlan->AddAblationZone(m_TempAblationZones.at(0));

      while (indices.size() > 0 &&
             (double)(indices.size() / size) >
               (m_ToleranceNonAblatedTumorSafetyMarginVolume / 100))
      {
        AblationUtils::AblationZone newAblationCenter =
          AblationUtils::SearchNextAblationCenter(indices,
                                                  currentPlan->GetSegmentationImage(),
                                                  m_AblationRadius,
                                                  m_MaxAblationRadius,
                                                  m_ImageDimension,
                                                  m_ImageSpacing);
        AblationUtils::CalculateAblationVolume(newAblationCenter.indexCenter,
                                               currentPlan->GetSegmentationImage(),
                                               newAblationCenter.radius,
                                               m_ImageSpacing,
                                               m_ImageDimension,
                                               m_TempAblationZones);
        AblationUtils::RemoveAblatedPixelsFromGivenVector(newAblationCenter.indexCenter,
                                                          indices,
                                                          currentPlan->GetSegmentationImage(),
                                                          newAblationCenter.radius,
                                                          m_ImageDimension,
                                                          m_ImageSpacing);
        currentPlan->AddAblationZone(newAblationCenter);
      }
    }
    //------------ End calculation models -------------------------------

    // Check if the radius of some ablation zones can be reduced
    for (int i = 0; i < currentPlan->GetNumberOfZones(); i++)
    {
      AblationUtils::AblationZone *zone = currentPlan->GetAblationZone(i);
      double currentRadius = AblationUtils::FindMinimalAblationRadius(zone->indexCenter,
                                                                      currentPlan->GetSegmentationImage(),
                                                                      zone->radius,
                                                                      m_MinAblationRadius,
                                                                      m_ImageDimension,
                                                                      m_ImageSpacing);
      // MITK_INFO << "Found minimal radius: " << currentRadius;
      (*zone).radius = currentRadius;
    }

    // Check if some zones can be removed (TODO: fix!)
    //currentPlan->DetectAndRemoveNotNeededVolumes();

    MITK_INFO << "Total number of ablation zones: " << currentPlan->GetNumberOfZones();
    AllFoundPlans.push_back(currentPlan);
    //AllFoundPlans[iteration] = currentPlan;
  } // End of for loop

  // TODO
  // Remove vectors that are not needed any more
  // this->CopyTemporaryAblationZoneDistribution(); (REMOVE METHOD!)
  // AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension); COPY IMAGE
  MITK_INFO << "Found " << AllFoundPlans.size() << " proposals";
  // Search for best final proposal:
  mitk::AblationPlan::Pointer finalProposal = AllFoundPlans.at(0);
  for (int i = 1; i < AllFoundPlans.size(); i++)
  {
    if (finalProposal->CompareTo(AllFoundPlans.at(i)) == 1)
    {
      finalProposal = AllFoundPlans.at(i);
      MITK_INFO << "Best proposal: " << i;
    }
  }

  //==================== Optimization of final proposal ==================================================

  // Check if ablation zones are too far outside the tumor, if yes move them towards the center
  for (int index = 0; index < finalProposal->GetNumberOfZones(); ++index)
  {
    double ratio = AblationUtils::CalculateRatioAblatedTissueOutsideTumorToAblatedTissueInsideTumor(
      finalProposal->GetAblationZone(index)->indexCenter,
      finalProposal->GetSegmentationImage(),
      finalProposal->GetAblationZone(index)->radius,
      finalProposal->GetImageDimension(),
      finalProposal->GetImageSpacing());
    MITK_WARN << "RATIO: " << ratio;
    if (ratio > 0.3)
    {
      AblationUtils::MoveCenterOfAblationZone(finalProposal->GetAblationZone(index)->indexCenter,
                                              finalProposal->GetSegmentationImage(),
                                              finalProposal->GetAblationZone(index)->radius,
                                              finalProposal->GetImageDimension(),
                                              finalProposal->GetImageSpacing());
    }
  }

  // Check, if ablation zones have a too short distance between each other, if yes they can be removed
  for (int index = 0; index < finalProposal->GetNumberOfZones(); ++index)
  {
    std::vector<int> indexToRemove;
    itk::Index<3> actualIndex = finalProposal->GetAblationZone(index)->indexCenter;
    for (int counter = 0; counter < finalProposal->GetNumberOfZones(); ++counter)
    {
      if (counter == index)
      {
        continue;
      }
      itk::Index<3> indexToProof = finalProposal->GetAblationZone(counter)->indexCenter;
      double distance = AblationUtils::CalculateScalarDistance(actualIndex, indexToProof, m_ImageSpacing);
      if (distance <= 0.5 * finalProposal->GetAblationZone(counter)->radius)
      {
        indexToRemove.push_back(counter);
      }
    }
    for (int position = indexToRemove.size() - 1; position >= 0; --position)
    {
      finalProposal->RemoveAblationZone(indexToRemove.at(position));
      MITK_DEBUG << "Removed Ablation zone at index position: " << indexToRemove.at(position);
      index = -1;
    }
  }



  for (int index = 0; index < finalProposal->GetNumberOfZones(); ++index)
  {
    AblationUtils::CalculateAblationVolume(finalProposal->GetAblationZone(index)->indexCenter,
                                           m_SegmentationImage,
                                           finalProposal->GetAblationZone(index)->radius,
                                           m_ImageSpacing,
                                           m_ImageDimension,
                                           m_TempAblationZones);
  }

  /*TODO: Adapt method...
  AblationUtils::DetectNotNeededAblationVolume(
    m_AblationZonesProcessed, m_AblationZones, m_SegmentationImage, m_ImageDimension, m_ImageSpacing);*/



  //============
  /* Todo: Adapt!
  if (m_Controls.toleranceNonAblatedTumorSafetyMarginVolumeSpinBox->value() == 0)
  {
    std::vector<itk::Index<3>> onlyTumorIndices =
      AblationUtils::FillVectorContainingIndicesOfTumorTissueOnly(m_SegmentationImage, m_ImageDimension);

    while (onlyTumorIndices.size() > 0)
    {
      AblationUtils::AblationZone newAblationCenter = AblationUtils::SearchNextAblationCenter(
        onlyTumorIndices, m_SegmentationImage, m_AblationRadius, m_MaxAblationRadius, m_ImageDimension, m_ImageSpacing);
      AblationUtils::MoveCenterOfAblationZone(
        newAblationCenter.indexCenter, m_SegmentationImage, newAblationCenter.radius, m_ImageDimension, m_ImageSpacing);
      AblationUtils::CalculateAblationVolume(
        newAblationCenter.indexCenter, m_SegmentationImage, newAblationCenter.radius, m_ImageSpacing, m_ImageDimension);
      AblationUtils::RemoveAblatedPixelsFromGivenVector(newAblationCenter.indexCenter,
                                                        onlyTumorIndices,
                                                        m_SegmentationImage,
                                                        newAblationCenter.radius,
                                                        m_ImageDimension,
                                                        m_ImageSpacing);
      m_AblationZonesProcessed.push_back(newAblationCenter);
      m_AblationZones.push_back(newAblationCenter);
    }
    AblationUtils::DetectNotNeededAblationVolume(
      m_AblationZonesProcessed, m_AblationZones, m_SegmentationImage, m_ImageDimension, m_ImageSpacing);
  }
  //============
  */
  m_AblationPlan = finalProposal;
}
