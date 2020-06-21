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
#include "mitkAblationPlanOptimizer.h"

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
  AllFoundPlans.resize(m_Iterations);

  // Start of for-loop (main iterative loop, each iteration = one proposal):
  MITK_INFO << " Creating " << m_Iterations << " proposal templates ...";
  std::vector<mitk::AblationPlan::Pointer> planTemplates;
  for(int i=0; i<m_Iterations; i++){
    mitk::AblationPlan::Pointer currentPlan = mitk::AblationPlan::New();
    currentPlan->SetSegmentationImage(m_SegmentationImage->Clone());
    currentPlan->SetImageDimension(m_ImageDimension);
    currentPlan->SetImageSpacing(m_ImageSpacing);
    planTemplates.push_back(currentPlan);
    if (i%25 == 0) {
      MITK_INFO << "..." << i;
    }
  }
  MITK_INFO << " Computing " << m_Iterations << " proposal in parallel (CPU) ...";
  #pragma omp parallel for
  for (int iteration = 1; iteration <= m_Iterations; ++iteration)
  {
    mitk::AblationPlan::Pointer currentPlan = planTemplates[iteration-1];
    std::vector<mitk::AblationZone> tempAblationZones;

    std::vector<itk::Index<3>> tumorTissueSafetyMarginIndices;
    for(itk::Index<3> i : m_TumorTissueSafetyMarginIndices){
      tumorTissueSafetyMarginIndices.push_back(i);
    }

    //if no manual starting point is set: find new starting point
    //if (!m_ManualAblationStartingPositionSet){
      QString position = AblationUtils::FindAblationStartingPosition(currentPlan->GetSegmentationImage(),
                                                                     tumorTissueSafetyMarginIndices,
                                                                     m_AblationRadius,
                                                                     m_MaxAblationRadius,
                                                                     m_TempAblationStartingPositionIndexCoordinates,
                                                                     m_TempAblationStartingPositionInWorldCoordinates,
                                                                     m_AblationRadius,
                                                                     m_ImageDimension,
                                                                     m_ImageSpacing);
      //m_Controls.ablationStartingPointLabel->setText(position);
      //MITK_INFO << "Found starting point: " << position.toLatin1().toStdString();
    //}


    MITK_INFO << "... new iteration";


    //------------ Random distribution model calculations: ---------------
    {
      double size = tumorTissueSafetyMarginIndices.size();
      std::vector<itk::Index<3>> indices = tumorTissueSafetyMarginIndices;

      // Add starting zone
      AblationUtils::CalculateAblationVolume(m_TempAblationStartingPositionIndexCoordinates,
                                             currentPlan->GetSegmentationImage(),
                                             m_AblationRadius,
                                             m_ImageSpacing,
                                             m_ImageDimension,
                                             tempAblationZones);
      AblationUtils::RemoveAblatedPixelsFromGivenVector(tempAblationZones.at(0).indexCenter,
                                                        indices,
                                                        currentPlan->GetSegmentationImage(),
                                                        tempAblationZones.at(0).radius,
                                                        m_ImageDimension,
                                                        m_ImageSpacing);
      currentPlan->AddAblationZone(tempAblationZones.at(0));
      // end add starting zone

      while (indices.size() > 0 &&
             (double)(indices.size() / size) >
               (m_ToleranceNonAblatedTumorSafetyMarginVolume / 100))
      {
        mitk::AblationZone newAblationCenter =
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
                                               tempAblationZones);
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

    /* Another try to improve algorithm: Check if the radius of some ablation zones can be reduced
    for (int i = 0; i < currentPlan->GetNumberOfZones(); i++)
    {
      mitk::AblationZone *zone = currentPlan->GetAblationZone(i);
      double currentRadius = AblationUtils::FindMinimalAblationRadius(zone->indexCenter,
                                                                      currentPlan->GetSegmentationImage(),
                                                                      zone->radius,
                                                                      m_MinAblationRadius,
                                                                      m_ImageDimension,
                                                                      m_ImageSpacing);
      // MITK_INFO << "Found minimal radius: " << currentRadius;
      (*zone).radius = currentRadius;
    } */

    //MITK_INFO << "Number of ablation zones before reduction: " << currentPlan->GetNumberOfZones();

    //Optimize this proposal
    mitk::AblationPlanOptimizer::Optimize(currentPlan,tempAblationZones);

    // Check if some zones can be removed
    //AblationUtils::DetectNotNeededAblationVolume(currentPlan,currentPlan->GetSegmentationImage(),currentPlan->GetImageDimension(),currentPlan->GetImageSpacing());

    //MITK_INFO << "Final number of ablation zones: " << currentPlan->GetNumberOfZones();
    AblationUtils::ComputeStatistics(currentPlan,m_TumorTissueSafetyMarginIndices);
    AllFoundPlans[iteration-1] = currentPlan;
  } // End of for loop

  // TODO
  // Remove vectors that are not needed any more
  // this->CopyTemporaryAblationZoneDistribution(); (REMOVE METHOD!)
  // AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension); COPY IMAGE
  MITK_INFO << "Found " << AllFoundPlans.size() << " proposals";
  // Search for best final proposal:
  mitk::AblationPlan::Pointer finalProposal = AllFoundPlans.at(0);
  //MITK_INFO << "Proposal 0: zones: " << AllFoundPlans.at(0)->GetNumberOfZones() <<"; overlap: " << AllFoundPlans.at(0)->GetStatistics().factorOverlappingAblationZones;
  for (int i = 1; i < AllFoundPlans.size(); i++)
  {
    //MITK_INFO << "Proposal " << i << ": zones: " << AllFoundPlans.at(i)->GetNumberOfZones() <<"; overlap: " << AllFoundPlans.at(i)->GetStatistics().factorOverlappingAblationZones;
    if (AllFoundPlans.at(i)->CompareTo(finalProposal) == -1)
    {
      finalProposal = AllFoundPlans.at(i);
      MITK_INFO << "Best proposal: " << i << " ("<<finalProposal->GetNumberOfZones()<<";"<<finalProposal->GetStatistics().factorNonAblatedVolume<<")";
    }
  }

  //==================== Optimization of final proposal ==================================================
  /*
  mitk::AblationPlanOptimizer::Optimize(finalProposal,tempAblationZones);
  //============

  if (m_ToleranceNonAblatedTumorSafetyMarginVolume == 0)
  {
    std::vector<itk::Index<3>> onlyTumorIndices =
      AblationUtils::FillVectorContainingIndicesOfTumorTissueOnly(finalProposal->GetSegmentationImage(), m_ImageDimension);

    while (onlyTumorIndices.size() > 0)
    {
      mitk::AblationZone newAblationCenter = AblationUtils::SearchNextAblationCenter(
        onlyTumorIndices, finalProposal->GetSegmentationImage(), m_AblationRadius, m_MaxAblationRadius, m_ImageDimension, m_ImageSpacing);
      AblationUtils::MoveCenterOfAblationZone(
        newAblationCenter.indexCenter, finalProposal->GetSegmentationImage(), newAblationCenter.radius, m_ImageDimension, m_ImageSpacing);
      AblationUtils::CalculateAblationVolume(
        newAblationCenter.indexCenter, finalProposal->GetSegmentationImage(), newAblationCenter.radius, m_ImageSpacing, m_ImageDimension);
      AblationUtils::RemoveAblatedPixelsFromGivenVector(newAblationCenter.indexCenter,
                                                        onlyTumorIndices,
                                                        m_SegmentationImage,
                                                        newAblationCenter.radius,
                                                        m_ImageDimension,
                                                        m_ImageSpacing);
      finalProposal->AddAblationZone(newAblationCenter);
    }
    AblationUtils::DetectNotNeededAblationVolume(
     finalProposal, finalProposal->GetSegmentationImage(), m_ImageDimension, m_ImageSpacing);
  }*/
  //============
  m_AblationPlan = finalProposal;
}
