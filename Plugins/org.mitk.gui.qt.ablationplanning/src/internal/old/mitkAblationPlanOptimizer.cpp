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

#include "mitkAblationPlanOptimizer.h"
#include "mitkAblationUtils.h"

mitk::AblationPlanOptimizer::AblationPlanOptimizer() {}

mitk::AblationPlanOptimizer::~AblationPlanOptimizer() {}

void mitk::AblationPlanOptimizer::RemoveNotNeededVolumes(mitk::AblationPlan::Pointer plan){
  AblationUtils::DetectNotNeededAblationVolume(plan,plan->GetSegmentationImage(),plan->GetImageDimension(),plan->GetImageSpacing());
}

void mitk::AblationPlanOptimizer::Optimize(mitk::AblationPlan::Pointer plan, std::vector<mitk::AblationZone> &tempAblationZones){
  //==================== Optimization of final proposal ==================================================
  AblationUtils::ResetSegmentationImage(plan->GetSegmentationImage(), plan->GetImageDimension());

  // Check if ablation zones are too far outside the tumor, if yes move them towards the center
  for (int index = 0; index < plan->GetNumberOfZones(); ++index)
  {
    double ratio = AblationUtils::CalculateRatioAblatedTissueOutsideTumorToAblatedTissueInsideTumor(
      plan->GetAblationZone(index)->indexCenter,
      plan->GetSegmentationImage(),
      plan->GetAblationZone(index)->radius,
      plan->GetImageDimension(),
      plan->GetImageSpacing());
    //MITK_WARN << "RATIO: " << ratio;
    if (ratio > 0.3)
    {
      AblationUtils::MoveCenterOfAblationZone(plan->GetAblationZone(index)->indexCenter,
                                              plan->GetSegmentationImage(),
                                              plan->GetAblationZone(index)->radius,
                                              plan->GetImageDimension(),
                                              plan->GetImageSpacing());
    }
  }


  // Check, if ablation zones have a too short distance between each other, if yes they can be removed
  for (int index = 0; index < plan->GetNumberOfZones(); ++index)
  {
    std::vector<int> indexToRemove;
    itk::Index<3> actualIndex = plan->GetAblationZone(index)->indexCenter;
    for (int counter = 0; counter < plan->GetNumberOfZones(); ++counter)
    {
      if (counter == index)
      {
        continue;
      }
      itk::Index<3> indexToProof = plan->GetAblationZone(counter)->indexCenter;
      double distance = AblationUtils::CalculateScalarDistance(actualIndex, indexToProof, plan->GetImageSpacing());
      if (distance <= 0.5 * plan->GetAblationZone(counter)->radius)
      {
        indexToRemove.push_back(counter);
      }
    }
    for (int position = indexToRemove.size() - 1; position >= 0; --position)
    {
      plan->RemoveAblationZone(indexToRemove.at(position));
      MITK_DEBUG << "Removed Ablation zone at index position: " << indexToRemove.at(position);
      index = -1;
    }
  }

  //Detect and delet zones that are not needed any more after the others where moved
  for (int index = 0; index < plan->GetNumberOfZones(); ++index)
  {
    AblationUtils::CalculateAblationVolume(plan->GetAblationZone(index)->indexCenter,
                                           plan->GetSegmentationImage(),
                                           plan->GetAblationZone(index)->radius,
                                           plan->GetImageSpacing(),
                                           plan->GetImageDimension(),
                                           tempAblationZones);
  }

  RemoveNotNeededVolumes(plan);





  //============
}
