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

const static short ABLATION_VALUE = 2;
const static short TUMOR_NOT_YET_ABLATED = 1;
const static short NO_TUMOR_ISSUE = 0;
const static short SAFETY_MARGIN = 256;
const static unsigned short BIT_OPERATION_ELIMINATE_TUMOR_SAFETY_MARGIN = 65278; // = 11111110 11111110

mitk::AblationPlanOptimizer::AblationPlanOptimizer() {}

mitk::AblationPlanOptimizer::~AblationPlanOptimizer() {}


//void mitk::AblationPlanOptimizer::SetNewAblationVolume(mitk::AblationPlan::Pointer plan,
//                                                       int zoneNumber,
//                                                       itk::Index<3> newCenter)
//{
//  mitk::Image::Pointer image = plan->GetSegmentationImage();
//  mitk::Vector3D imageSpacing = plan->GetImageSpacing();
//  mitk::Vector3D imageDimension = plan->GetImageDimension();
//  double radius = plan->GetAblationZone(zoneNumber)->radius;
//  if (image.IsNotNull())
//  {
//    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
//    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
//    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);
//
//    unsigned int upperX;
//    unsigned int lowerX;
//    unsigned int upperY;
//    unsigned int lowerY;
//    unsigned int upperZ;
//    unsigned int lowerZ;
//
//    AblationUtils::CalculateUpperLowerXYZ(upperX,
//                                          lowerX,
//                                          upperY,
//                                          lowerY,
//                                          upperZ,
//                                          lowerZ,
//                                          pixelDirectionX,
//                                          pixelDirectionY,
//                                          pixelDirectionZ,
//                                          newCenter,
//                                          imageDimension);
//
//    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
//    itk::Index<3> actualIndex;
//    int pixelValue{0};
//    for (actualIndex[2] = lowerZ - 1; actualIndex[2] <= upperZ + 1; actualIndex[2] += 1)
//    {
//      for (actualIndex[1] = lowerY - 1; actualIndex[1] <= upperY + 1; actualIndex[1] += 1)
//      {
//        for (actualIndex[0] = lowerX - 1; actualIndex[0] <= upperX + 1; actualIndex[0] += 1)
//        {
//          if (radius <= AblationUtils::CalculateScalarDistance(newCenter, actualIndex, imageSpacing) &&
//              radius > AblationUtils::CalculateScalarDistance(
//                          plan->GetAblationZone(zoneNumber)->indexCenter, actualIndex, imageSpacing))
//          {
//            pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
//            imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue - ABLATION_VALUE);
//          }
//        }
//      }
//    }
//  }
//}

//bool mitk::AblationPlanOptimizer::CheckNewCenterForNonAblatedVolume(mitk::AblationPlan::Pointer plan,
//                                                                    int zoneNumber,
//                                                                    itk::Index<3> newCenter)
//{
//  mitk::Image::Pointer image = plan->GetSegmentationImage();
//  mitk::Vector3D imageSpacing = plan->GetImageSpacing();
//  mitk::Vector3D imageDimension = plan->GetImageDimension();
//  double radius = plan->GetAblationZone(zoneNumber)->radius;
//  if (image.IsNotNull())
//  {
//    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
//    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
//    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);
//
//    unsigned int upperX;
//    unsigned int lowerX;
//    unsigned int upperY;
//    unsigned int lowerY;
//    unsigned int upperZ;
//    unsigned int lowerZ;
//
//    AblationUtils::CalculateUpperLowerXYZ(upperX,
//                                          lowerX,
//                                          upperY,
//                                          lowerY,
//                                          upperZ,
//                                          lowerZ,
//                                          pixelDirectionX,
//                                          pixelDirectionY,
//                                          pixelDirectionZ,
//                                          newCenter,
//                                          imageDimension);
//
//    mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
//    itk::Index<3> actualIndex;
//    int pixelValue{0};
//    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
//    {
//      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
//      {
//        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
//        {
//          if (radius <= AblationUtils::CalculateScalarDistance(newCenter, actualIndex, imageSpacing) &&
//              radius >= AblationUtils::CalculateScalarDistance(
//                          plan->GetAblationZone(zoneNumber)->indexCenter, actualIndex, imageSpacing))
//          {
//            pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
//            if (pixelValue - ABLATION_VALUE == TUMOR_NOT_YET_ABLATED || pixelValue - ABLATION_VALUE == SAFETY_MARGIN)
//            {
//              return false;
//            }
//          }
//        }
//      }
//    }
//    return true;
//  }
//  return false;
//}

void mitk::AblationPlanOptimizer::Optimize(mitk::AblationPlan::Pointer plan,
                                           std::vector<mitk::AblationZone> &tempAblationZones)
{
  //==================== Optimization of final proposal ==================================================
  AblationUtils::ResetSegmentationImage(plan->GetSegmentationImage(), plan->GetImageDimension());

  // Check if ablation zones are too far outside the tumor, if yes move them towards the center, OBSOLETE WITH THE NEW CHANGES
  //MoveZonesTowardsCenter(plan, tempAblationZones);

  // Check, if ablation zones have a too short distance between each other, if yes they can be removed, OBSOLETE WITH THE NEW CHANGES
  //CheckIfDistanceBeetweenZonesTooSmall(plan, tempAblationZones);

  // Detect and delet zones that are not needed any more after the others where moved
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

void mitk::AblationPlanOptimizer::RemoveNotNeededVolumes(mitk::AblationPlan::Pointer plan)
{
  AblationUtils::DetectNotNeededAblationVolume(
    plan, plan->GetSegmentationImage(), plan->GetImageDimension(), plan->GetImageSpacing());
}

void mitk::AblationPlanOptimizer::MoveZonesTowardsCenter(mitk::AblationPlan::Pointer plan,
                                                         std::vector<mitk::AblationZone> &tempAblationZones)
{
  for (int index = 0; index < plan->GetNumberOfZones(); ++index)
  {
    double ratio = AblationUtils::CalculateRatioAblatedTissueOutsideTumorToAblatedTissueInsideTumor(
      plan->GetAblationZone(index)->indexCenter,
      plan->GetSegmentationImage(),
      plan->GetAblationZone(index)->radius,
      plan->GetImageDimension(),
      plan->GetImageSpacing());
    // MITK_WARN << "RATIO: " << ratio;
    if (ratio > 0.3)
    {
      AblationUtils::MoveCenterOfAblationZone(plan->GetAblationZone(index)->indexCenter,
                                              plan->GetSegmentationImage(),
                                              plan->GetAblationZone(index)->radius,
                                              plan->GetImageDimension(),
                                              plan->GetImageSpacing());
    }
  }
}

void mitk::AblationPlanOptimizer::CheckIfDistanceBeetweenZonesTooSmall(
  mitk::AblationPlan::Pointer plan, std::vector<mitk::AblationZone> &tempAblationZones)
{
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
}