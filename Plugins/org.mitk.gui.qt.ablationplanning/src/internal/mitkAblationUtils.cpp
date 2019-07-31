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

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "mitkAblationUtils.h"

// Qt
#include <QMessageBox>

// mitk image
#include "mitkProperties.h"
#include <mitkImage.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>

#include <mitkSurface.h>
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>

#include <cmath>

const static short ABLATION_VALUE = 2;
const static short TUMOR_NOT_YET_ABLATED = 1;
const static short NO_TUMOR_ISSUE = 0;
const static short SAFETY_MARGIN = 256;
const static unsigned short BIT_OPERATION_ELIMINATE_TUMOR_SAFETY_MARGIN = 65278; // = 11111110 11111110

//=====================Konstruktor/Destruktor===================================
AblationUtils::AblationUtils() {}

AblationUtils::~AblationUtils() {}

void AblationUtils::FillVectorContainingIndicesOfTumorTissueSafetyMargin(
  mitk::Image::Pointer image,
  mitk::Vector3D &imageDimension,
  std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices)
{
  if (image.IsNotNull())
  {
    MITK_DEBUG << "Detecting the tumor tissue and safety margin is in progress...";
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          if (imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED ||
              imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN)
          {
            tumorTissueSafetyMarginIndices.push_back(actualIndex);
          }
        }
      }
    }
  }
}

std::vector<itk::Index<3>> AblationUtils::FillVectorContainingIndicesOfTumorTissueOnly(mitk::Image::Pointer image,
                                                                                       mitk::Vector3D &imageDimension)
{
  std::vector<itk::Index<3>> onlyTumorIndices;
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          if (imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED)
          {
            onlyTumorIndices.push_back(actualIndex);
          }
        }
      }
    }
  }
  return onlyTumorIndices;
}

QString AblationUtils::FindAblationStartingPosition(mitk::Image::Pointer image,
                                                    std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices,
                                                    double &ablationRadius,
                                                    itk::Index<3> &tempAblationStartingPositionIndexCoordinates,
                                                    mitk::Point3D &tempAblationStartingPositionInWorldCoordinates,
                                                    mitk::Vector3D &imageDimension,
                                                    mitk::Vector3D &imageSpacing)
{
  if (image.IsNotNull())
  {
    MITK_DEBUG << "Finding a random ablation starting position...";

    bool positionFound = false;
    int iteration = 1;

    while (!positionFound && iteration < 20)
    {
      int randomIndex1 = rand() % tumorTissueSafetyMarginIndices.size();
      int randomIndex2 = rand() % tumorTissueSafetyMarginIndices.size();
      int randomIndex3 = rand() % tumorTissueSafetyMarginIndices.size();
      int randomIndex4 = rand() % tumorTissueSafetyMarginIndices.size();
      int randomIndex5 = rand() % tumorTissueSafetyMarginIndices.size();

      itk::Index<3> startingPosition1 = tumorTissueSafetyMarginIndices.at(randomIndex1);
      itk::Index<3> startingPosition2 = tumorTissueSafetyMarginIndices.at(randomIndex2);
      itk::Index<3> startingPosition3 = tumorTissueSafetyMarginIndices.at(randomIndex3);
      itk::Index<3> startingPosition4 = tumorTissueSafetyMarginIndices.at(randomIndex4);
      itk::Index<3> startingPosition5 = tumorTissueSafetyMarginIndices.at(randomIndex5);

      std::vector<itk::Index<3>> startingPositions;
      startingPositions.push_back(startingPosition1);
      startingPositions.push_back(startingPosition2);
      startingPositions.push_back(startingPosition3);
      startingPositions.push_back(startingPosition4);
      startingPositions.push_back(startingPosition5);

      double radius1 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition1, image, imageSpacing, imageDimension);
      double radius2 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition2, image, imageSpacing, imageDimension);
      double radius3 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition3, image, imageSpacing, imageDimension);
      double radius4 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition4, image, imageSpacing, imageDimension);
      double radius5 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition5, image, imageSpacing, imageDimension);
      std::vector<double> radiusVector;
      std::vector<double>::iterator result;
      radiusVector.push_back(radius1);
      radiusVector.push_back(radius2);
      radiusVector.push_back(radius3);
      radiusVector.push_back(radius4);
      radiusVector.push_back(radius5);
      result = std::max_element(radiusVector.begin(), radiusVector.end());
      int index = std::distance(radiusVector.begin(), result);

      if (radiusVector.at(index) >= ablationRadius - 1)
      {
        positionFound = true;
      }
      else if (iteration > 10 && radiusVector.at(index) < ablationRadius - 1 &&
               radiusVector.at(index) > ablationRadius / 2)
      {
        positionFound = true;
      }
      else
      {
        ++iteration;
      }
      MITK_DEBUG << "Iteration: " << iteration;
      tempAblationStartingPositionIndexCoordinates = startingPositions.at(index);
    }

    // Calculate the index coordinates of the starting position:
    image->GetGeometry()->IndexToWorld(tempAblationStartingPositionIndexCoordinates,
                                       tempAblationStartingPositionInWorldCoordinates);

    double x = tempAblationStartingPositionInWorldCoordinates[0];
    double y = tempAblationStartingPositionInWorldCoordinates[1];
    double z = tempAblationStartingPositionInWorldCoordinates[2];
    QString text = QString("Set Ablation Startingposition to: %1 | %2 | %3").arg(x).arg(y).arg(z);

    MITK_DEBUG << "Set Ablation Startingposition to: " << tempAblationStartingPositionInWorldCoordinates;
    MITK_DEBUG << "Startingposition in Index: " << tempAblationStartingPositionIndexCoordinates;
    MITK_DEBUG << "Spacing: " << image->GetGeometry()->GetSpacing();
    // Get number of voxels in the three dimensions:
    MITK_DEBUG << "Dimension: " << imageDimension[0] << " " << imageDimension[1] << " " << imageDimension[2];
    return text;
  }
  return "";
}

double AblationUtils::CalculateScalarDistance(itk::Index<3> &point1,
                                              itk::Index<3> &point2,
                                              mitk::Vector3D &imageSpacing)
{
  double x = (point1[0] - point2[0]) * imageSpacing[0];
  double y = (point1[1] - point2[1]) * imageSpacing[1];
  double z = (point1[2] - point2[2]) * imageSpacing[2];

  return sqrt(x * x + y * y + z * z);
}

void AblationUtils::CalculateAblationVolume(itk::Index<3> &center,
                                            mitk::Image::Pointer image,
                                            double &radius,
                                            mitk::Vector3D &imageSpacing,
                                            mitk::Vector3D &imageDimension,
                                            std::vector<itk::Index<3>> &tempAblationZoneCenters)
{
  MITK_DEBUG << "Calculate ablation volume for index: " << center;
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           center,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(center, actualIndex, imageSpacing))
          {
            unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex) + ABLATION_VALUE;

            imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
          }
        }
      }
    }
    tempAblationZoneCenters.push_back(center);
  }
}

void AblationUtils::CalculateAblationVolume(itk::Index<3> &center,
                                            mitk::Image::Pointer image,
                                            double &radius,
                                            mitk::Vector3D &imageSpacing,
                                            mitk::Vector3D &imageDimension)
{
  MITK_DEBUG << "Calculate ablation volume for index: " << center;
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           center,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(center, actualIndex, imageSpacing))
          {
            unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex) + ABLATION_VALUE;

            imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
          }
        }
      }
    }
  }
}

bool AblationUtils::CheckVolumeForNonAblatedTissue(itk::Index<3> &centerOfVolume,
                                                   mitk::Image::Pointer image,
                                                   double &radius,
                                                   mitk::Vector3D &imageSpacing,
                                                   mitk::Vector3D &imageDimension)
{
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           centerOfVolume,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(centerOfVolume, actualIndex, imageSpacing))
          {
            if (imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED ||
                imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN)
            {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

bool AblationUtils::CheckIfVolumeOfGivenRadiusIsTotallyInsideTumorTissueAndSafetyMargin(double &radius,
                                                                                        itk::Index<3> &centerOfVolume,
                                                                                        mitk::Image::Pointer image,
                                                                                        mitk::Vector3D &imageSpacing,
                                                                                        mitk::Vector3D &imageDimension)
{
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           centerOfVolume,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(centerOfVolume, actualIndex, imageSpacing))
          {
            if (imagePixelWriter.GetPixelByIndex(actualIndex) != TUMOR_NOT_YET_ABLATED &&
                imagePixelWriter.GetPixelByIndex(actualIndex) != SAFETY_MARGIN)
            {
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}

double AblationUtils::CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(itk::Index<3> &point,
                                                                         mitk::Image::Pointer image,
                                                                         mitk::Vector3D &imageSpacing,
                                                                         mitk::Vector3D &imageDimension)
{
  double radius = 1.0;
  while (CheckIfVolumeOfGivenRadiusIsTotallyInsideTumorTissueAndSafetyMargin(
    radius, point, image, imageSpacing, imageDimension))
  {
    radius += 1;
  }
  ++radius;
  // MITK_INFO << "Calculated max radius for given point: " << radius;
  return radius;
}

bool AblationUtils::CheckImageForNonAblatedTissue(mitk::Image::Pointer image, mitk::Vector3D &imageDimension)
{
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          if (imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED ||
              imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN)
          {
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool AblationUtils::CheckForNonAblatedTumorTissueWithoutSafetyMargin(std::vector<itk::Index<3>> &indices,
                                                                     mitk::Image::Pointer image,
                                                                     mitk::Vector3D &imageDimension)
{
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          if (imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED)
          {
            return true;
          }
        }
      }
    }
  }
  return false;
}

void AblationUtils::ProcessDirectNeighbourAblationZones(itk::Index<3> &center,
                                                        mitk::Image::Pointer image,
                                                        mitk::Vector3D &imageSpacing,
                                                        mitk::Vector3D &imageDimension,
                                                        double &ablationRadius,
                                                        std::vector<itk::Index<3>> &tempAblationZoneCentersProcessed,
                                                        std::vector<itk::Index<3>> &tempAblationZoneCenters)
{
  MITK_INFO << "Process direct neighbour ablation zones for index: " << center;
  std::vector<itk::Index<3>> indices =
    CalculateIndicesOfDirectNeighbourAblationZones(center, ablationRadius, imageSpacing, imageDimension);

  for (std::vector<itk::Index<3>>::iterator it = indices.begin(); it != indices.end(); ++it)
  {
    if (CheckVolumeForNonAblatedTissue(*it, image, ablationRadius, imageSpacing, imageDimension))
    {
      CalculateAblationVolume(*it, image, ablationRadius, imageSpacing, imageDimension, tempAblationZoneCenters);
    }
  }

  // Now, all 12 direct neighbour ablation zones are processed. So add the
  // index of the given center to the processed ablation centers:
  tempAblationZoneCentersProcessed.push_back(center);
}

void AblationUtils::CalculateUpperLowerXYZ(unsigned int &upperX,
                                           unsigned int &lowerX,
                                           unsigned int &upperY,
                                           unsigned int &lowerY,
                                           unsigned int &upperZ,
                                           unsigned int &lowerZ,
                                           unsigned int &pixelDirectionX,
                                           unsigned int &pixelDirectionY,
                                           unsigned int &pixelDirectionZ,
                                           itk::Index<3> &center,
                                           mitk::Vector3D &imageDimension)
{
  // Calculate upperX --> means vector in direction [1,0,0]:
  if (center[0] + pixelDirectionX >= imageDimension[0])
  {
    upperX = imageDimension[0] - 1;
  }
  else
  {
    upperX = center[0] + pixelDirectionX;
  }
  // Calculate lowerX --> means vector in direction [-1,0,0]:
  if (center[0] - pixelDirectionX < 0)
  {
    lowerX = 0;
  }
  else
  {
    lowerX = center[0] - pixelDirectionX;
  }

  // Calculate upperY --> means vector in direction [0,1,0]:
  if (center[1] + pixelDirectionY >= imageDimension[1])
  {
    upperY = imageDimension[1] - 1;
  }
  else
  {
    upperY = center[1] + pixelDirectionY;
  }
  // Calculate lowerY --> means vector in direction [0,-1,0]:
  if (center[1] - pixelDirectionY < 0)
  {
    lowerY = 0;
  }
  else
  {
    lowerY = center[1] - pixelDirectionY;
  }

  // Calculate upperZ --> means vector in direction [0,0,1]:
  if (center[2] + pixelDirectionZ >= imageDimension[2])
  {
    upperZ = imageDimension[2] - 1;
  }
  else
  {
    upperZ = center[2] + pixelDirectionZ;
  }

  // Calculate lowerZ --> means vector in direction [0,0,-1]:
  if (center[2] - pixelDirectionZ < 0)
  {
    lowerZ = 0;
  }
  else
  {
    lowerZ = center[2] - pixelDirectionZ;
  }
}

void AblationUtils::CalculateDistancesOfTumorBoundariesFromCenter(double &distanceLowerX,
                                                                  double &distanceUpperX,
                                                                  double &distanceLowerY,
                                                                  double &distanceUpperY,
                                                                  double &distanceLowerZ,
                                                                  double &distanceUpperZ,
                                                                  itk::Index<3> &center,
                                                                  mitk::Image::Pointer image,
                                                                  mitk::Vector3D &imageDimension,
                                                                  mitk::Vector3D &imageSpacing)
{
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> newIndex = center;
    // Calculate distance from center to upper boundary in x direction:
    for (int count = center[0] + 1; count < imageDimension[0]; ++count)
    {
      newIndex[0] = count;
      if ((imagePixelWriter.GetPixelByIndex(newIndex) & TUMOR_NOT_YET_ABLATED) == TUMOR_NOT_YET_ABLATED ||
          (imagePixelWriter.GetPixelByIndex(newIndex) & SAFETY_MARGIN) == SAFETY_MARGIN)
      {
        continue;
      }
      else
      {
        distanceUpperX = CalculateScalarDistance(center, newIndex, imageSpacing);
        break;
      }
    }

    // Calculate distance from center to lower boundary in x direction:
    newIndex = center;
    for (int count = center[0] - 1; count >= 0; --count)
    {
      newIndex[0] = count;
      if ((imagePixelWriter.GetPixelByIndex(newIndex) & TUMOR_NOT_YET_ABLATED) == TUMOR_NOT_YET_ABLATED ||
          (imagePixelWriter.GetPixelByIndex(newIndex) & SAFETY_MARGIN) == SAFETY_MARGIN)
      {
        continue;
      }
      else
      {
        distanceLowerX = CalculateScalarDistance(center, newIndex, imageSpacing);
        break;
      }
    }

    // Calculate distance from center to upper boundary in y direction:
    newIndex = center;
    for (int count = center[1] + 1; count < imageDimension[1]; ++count)
    {
      newIndex[1] = count;
      if ((imagePixelWriter.GetPixelByIndex(newIndex) & TUMOR_NOT_YET_ABLATED) == TUMOR_NOT_YET_ABLATED ||
          (imagePixelWriter.GetPixelByIndex(newIndex) & SAFETY_MARGIN) == SAFETY_MARGIN)
      {
        continue;
      }
      else
      {
        distanceUpperY = CalculateScalarDistance(center, newIndex, imageSpacing);
        break;
      }
    }

    // Calculate distance from center to lower boundary in y direction:
    newIndex = center;
    for (int count = center[1] - 1; count >= 0; --count)
    {
      newIndex[1] = count;
      if ((imagePixelWriter.GetPixelByIndex(newIndex) & TUMOR_NOT_YET_ABLATED) == TUMOR_NOT_YET_ABLATED ||
          (imagePixelWriter.GetPixelByIndex(newIndex) & SAFETY_MARGIN) == SAFETY_MARGIN)
      {
        continue;
      }
      else
      {
        distanceLowerY = CalculateScalarDistance(center, newIndex, imageSpacing);
        break;
      }
    }
    // Calculate distance from center to upper boundary in z direction:
    newIndex = center;
    for (int count = center[2] + 1; count < imageDimension[2]; ++count)
    {
      newIndex[2] = count;
      if ((imagePixelWriter.GetPixelByIndex(newIndex) & TUMOR_NOT_YET_ABLATED) == TUMOR_NOT_YET_ABLATED ||
          (imagePixelWriter.GetPixelByIndex(newIndex) & SAFETY_MARGIN) == SAFETY_MARGIN)
      {
        continue;
      }
      else
      {
        distanceUpperZ = CalculateScalarDistance(center, newIndex, imageSpacing);
        break;
      }
    }

    // Calculate distance from center to lower boundary in z direction:
    newIndex = center;
    for (int count = center[2] - 1; count >= 0; --count)
    {
      newIndex[2] = count;
      if ((imagePixelWriter.GetPixelByIndex(newIndex) & TUMOR_NOT_YET_ABLATED) == TUMOR_NOT_YET_ABLATED ||
          (imagePixelWriter.GetPixelByIndex(newIndex) & SAFETY_MARGIN) == SAFETY_MARGIN)
      {
        continue;
      }
      else
      {
        distanceLowerZ = CalculateScalarDistance(center, newIndex, imageSpacing);
        break;
      }
    }
  }
}

std::vector<itk::Index<3>> AblationUtils::CalculateIndicesOfDirectNeighbourAblationZones(itk::Index<3> &center,
                                                                                         double &radius,
                                                                                         mitk::Vector3D &imageSpacing,
                                                                                         mitk::Vector3D &imageDimension)
{
  MITK_INFO << "Calculate indices of direct neighbour ablation zones...";
  std::vector<itk::Index<3>> directNeighbourAblationZones;
  unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
  unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
  unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

  unsigned int upperX;
  unsigned int lowerX;
  unsigned int upperY;
  unsigned int lowerY;
  unsigned int upperZ;
  unsigned int lowerZ;

  itk::Index<3> newIndex;

  CalculateUpperLowerXYZ(upperX,
                         lowerX,
                         upperY,
                         lowerY,
                         upperZ,
                         lowerZ,
                         pixelDirectionX,
                         pixelDirectionY,
                         pixelDirectionZ,
                         center,
                         imageDimension);

  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [1, 1, 0]:
  newIndex[0] = upperX;
  newIndex[1] = upperY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [1, 1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [-1, 1, 0]:
  newIndex[0] = lowerX;
  newIndex[1] = upperY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [-1, 1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [1, -1, 0]:
  newIndex[0] = upperX;
  newIndex[1] = lowerY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [1, -1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [-1, -1, 0]:
  newIndex[0] = lowerX;
  newIndex[1] = lowerY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [-1, -1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [1, 0, 1]:
  newIndex[0] = upperX;
  newIndex[1] = center[1];
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [1, 0, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [1, 0, -1]:
  newIndex[0] = upperX;
  newIndex[1] = center[1];
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [1, 0, -1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [-1, 0, 1]:
  newIndex[0] = lowerX;
  newIndex[1] = center[1];
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [-1, 0, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [-1, 0, -1]:
  newIndex[0] = lowerX;
  newIndex[1] = center[1];
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [-1, 0, -1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [0, 1, 1]:
  newIndex[0] = center[0];
  newIndex[1] = upperY;
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [0, 1, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [0, 1, -1]:
  newIndex[0] = center[0];
  newIndex[1] = upperY;
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [0, 1, -1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [0, -1, 1]:
  newIndex[0] = center[0];
  newIndex[1] = lowerY;
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [0, -1, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  // Calculate position in vector direction [0, -1, -1]:
  newIndex[0] = center[0];
  newIndex[1] = lowerY;
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  // MITK_INFO << "Index for [0, -1, -1] --> " << newIndex;

  return directNeighbourAblationZones;
}

bool AblationUtils::IsAblationZoneAlreadyProcessed(itk::Index<3> &center,
                                                   std::vector<itk::Index<3>> &tempAblationZoneCentersProcessed)
{
  for (std::vector<itk::Index<3>>::iterator it = tempAblationZoneCentersProcessed.begin();
       it != tempAblationZoneCentersProcessed.end();
       ++it)
  {
    if (center == (*it))
    {
      return true;
    }
  }
  return false;
}

void AblationUtils::DetectNotNeededAblationVolume(std::vector<itk::Index<3>> &tempAblationZoneCentersProcessed,
                                                  std::vector<itk::Index<3>> &tempAblationZoneCenters,
                                                  mitk::Image::Pointer image,
                                                  double &radius,
                                                  mitk::Vector3D &imageDimension,
                                                  mitk::Vector3D &imageSpacing)
{
  std::vector<int> indicesRemoved;
  for (int index = 0; index < tempAblationZoneCentersProcessed.size(); ++index)
  {
    if (!CheckIfAblationVolumeIsNeeded(
          tempAblationZoneCentersProcessed.at(index), image, radius, imageDimension, imageSpacing))
    {
      RemoveAblationVolume(tempAblationZoneCentersProcessed.at(index), image, radius, imageDimension, imageSpacing);
      indicesRemoved.push_back(index);
    }
  }
  for (int index = indicesRemoved.size() - 1; index >= 0; --index)
  {
    std::vector<itk::Index<3>>::iterator it = tempAblationZoneCentersProcessed.begin();
    tempAblationZoneCentersProcessed.erase(it + indicesRemoved.at(index));
    std::vector<itk::Index<3>>::iterator it2 = tempAblationZoneCenters.begin();
    tempAblationZoneCenters.erase(it2 + indicesRemoved.at(index));
    MITK_INFO << "Removed Ablation zone at index position: " << indicesRemoved.at(index);
  }
}

bool AblationUtils::CheckIfAblationVolumeIsNeeded(itk::Index<3> &center,
                                                  mitk::Image::Pointer image,
                                                  double &radius,
                                                  mitk::Vector3D &imageDimension,
                                                  mitk::Vector3D &imageSpacing)
{
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           center,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(center, actualIndex, imageSpacing))
          {
            if (imagePixelWriter.GetPixelByIndex(actualIndex) - ABLATION_VALUE == TUMOR_NOT_YET_ABLATED ||
                imagePixelWriter.GetPixelByIndex(actualIndex) - ABLATION_VALUE == SAFETY_MARGIN)
            {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

double AblationUtils::FindMinimalAblationRadius(itk::Index<3> &center,
                                                mitk::Image::Pointer image,
                                                double &maxRadius,
                                                double &minRadius,
                                                mitk::Vector3D &imageDimension,
                                                mitk::Vector3D &imageSpacing)
{
  double currentRadius = maxRadius;
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(maxRadius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(maxRadius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(maxRadius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           center,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;

    while (currentRadius >= minRadius)
    {
      for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
      {
        for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
        {
          for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
          {
            if (maxRadius >= CalculateScalarDistance(center, actualIndex, imageSpacing) &&
                currentRadius <= CalculateScalarDistance(center, actualIndex, imageSpacing))
            {
              if (imagePixelWriter.GetPixelByIndex(actualIndex) - ABLATION_VALUE == TUMOR_NOT_YET_ABLATED ||
                  imagePixelWriter.GetPixelByIndex(actualIndex) - ABLATION_VALUE == SAFETY_MARGIN)
              {
                break;
              }
            }
          }
        }
      }
      currentRadius--;
    }
  }
  return currentRadius;
}

void AblationUtils::RemoveAblationVolume(itk::Index<3> &center,
                                         mitk::Image::Pointer image,
                                         double &radius,
                                         mitk::Vector3D &imageDimension,
                                         mitk::Vector3D &imageSpacing)
{
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           center,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(center, actualIndex, imageSpacing))
          {
            unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
            pixelValue -= ABLATION_VALUE;
            imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
          }
        }
      }
    }
  }
}

void AblationUtils::RemoveAblatedPixelsFromGivenVector(itk::Index<3> &center,
                                                       std::vector<itk::Index<3>> &tumorSafetyMarginPixels,
                                                       mitk::Image::Pointer image,
                                                       double &radius,
                                                       mitk::Vector3D &imageDimension,
                                                       mitk::Vector3D &imageSpacing)
{
  MITK_DEBUG << "Removing ablated pixels...";
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           center,
                           imageDimension);

    itk::Index<3> actualIndex;
    MITK_DEBUG << "Size of tumorPixels before: " << tumorSafetyMarginPixels.size();

    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(center, actualIndex, imageSpacing))
          {
            for (std::vector<itk::Index<3>>::iterator it = tumorSafetyMarginPixels.begin();
                 it != tumorSafetyMarginPixels.end();
                 ++it)
            {
              if ((*it) == actualIndex)
              {
                tumorSafetyMarginPixels.erase(it);
                break;
              }
            }
          }
        }
      }
    }
    MITK_DEBUG << "Size of tumorPixels after: " << tumorSafetyMarginPixels.size();
  }
}

itk::Index<3> AblationUtils::SearchNextAblationCenter(std::vector<itk::Index<3>> &tumorSafetyMarginPixels,
                                                      mitk::Image::Pointer image,
                                                      double &radius,
                                                      mitk::Vector3D &imageDimension,
                                                      mitk::Vector3D &imageSpacing)
{
  if (image.IsNotNull())
  {
    MITK_DEBUG << "Searching the next ablation center...";

    int iteration = 1;

    while (iteration < 20)
    {
      int randomIndex1 = rand() % tumorSafetyMarginPixels.size();
      int randomIndex2 = rand() % tumorSafetyMarginPixels.size();
      int randomIndex3 = rand() % tumorSafetyMarginPixels.size();
      int randomIndex4 = rand() % tumorSafetyMarginPixels.size();
      int randomIndex5 = rand() % tumorSafetyMarginPixels.size();

      itk::Index<3> position1 = tumorSafetyMarginPixels.at(randomIndex1);
      itk::Index<3> position2 = tumorSafetyMarginPixels.at(randomIndex2);
      itk::Index<3> position3 = tumorSafetyMarginPixels.at(randomIndex3);
      itk::Index<3> position4 = tumorSafetyMarginPixels.at(randomIndex4);
      itk::Index<3> position5 = tumorSafetyMarginPixels.at(randomIndex5);

      std::vector<itk::Index<3>> positions;
      positions.push_back(position1);
      positions.push_back(position2);
      positions.push_back(position3);
      positions.push_back(position4);
      positions.push_back(position5);

      double radius1 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(position1, image, imageSpacing, imageDimension);
      double radius2 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(position2, image, imageSpacing, imageDimension);
      double radius3 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(position3, image, imageSpacing, imageDimension);
      double radius4 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(position4, image, imageSpacing, imageDimension);
      double radius5 =
        CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(position5, image, imageSpacing, imageDimension);
      std::vector<double> radiusVector;
      std::vector<double>::iterator result;
      radiusVector.push_back(radius1);
      radiusVector.push_back(radius2);
      radiusVector.push_back(radius3);
      radiusVector.push_back(radius4);
      radiusVector.push_back(radius5);
      result = std::max_element(radiusVector.begin(), radiusVector.end());
      int index = std::distance(radiusVector.begin(), result);

      if (radiusVector.at(index) >= radius - 4)
      {
        return positions.at(index);
      }
      else
      {
        ++iteration;
      }
      /*if( iteration > 10)
      {
        itk::Index<3> position = positions.at(index);
        double radiusDifference = m_AblationRadius - radiusVector.at(index) + 1;
        std::vector<itk::Index<3>> neighbourPositions =
          this->CalculateIndicesOfDirectNeighbourAblationZones(position, radiusDifference);
        for( int count = 0; count < neighbourPositions.size(); ++count )
        {
          double tempRadius =
            this->CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(
                                            neighbourPositions.at(count));
          if (tempRadius >= m_AblationRadius - 4)
          {
            return neighbourPositions.at(count);
          }
        }
        for (int count = 0; count < neighbourPositions.size(); ++count)
        {
          double tempRadius =
            this->CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(
              neighbourPositions.at(count));
          if (tempRadius >= radiusVector.at(index))
          {
            return neighbourPositions.at(count);
          }
        }
      }*/
      if (iteration == 20)
      {
        return positions.at(index);
      }
    }
  }
  return itk::Index<3>();
}

void AblationUtils::ResetSegmentationImage(mitk::Image::Pointer image, mitk::Vector3D &imageDimension)
{
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
          pixelValue &= (SAFETY_MARGIN + TUMOR_NOT_YET_ABLATED);
          imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
        }
      }
    }
  }
}

void AblationUtils::ResetSafetyMargin(mitk::Image::Pointer image, mitk::Vector3D &imageDimension)
{
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
          pixelValue &= (SAFETY_MARGIN - 1);
          imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
        }
      }
    }
  }
}

bool AblationUtils::CheckAllVonNeumannNeighbourPixelsAreTumorTissue(itk::Index<3> &pixel,
                                                                    mitk::Image::Pointer image,
                                                                    mitk::Vector3D &imageDimension)
{
  if (image.IsNotNull())
  {
    // MITK_INFO << "CheckAllVonNeumannNeighbourPixels... " << pixel[0] << " " << pixel[1] << " " << pixel[2];
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    if ((imagePixelWriter.GetPixelByIndex(pixel) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED)
    {
      // If true --> actual pixel is no tumor tissue, so skip this pixel:
      return true;
    }

    unsigned int xLower = pixel[0];
    if (xLower != 0)
      --xLower;

    unsigned int xUpper = pixel[0];
    if (xUpper != imageDimension[0] - 1)
      ++xUpper;

    unsigned int yLower = pixel[1];
    if (yLower != 0)
      --yLower;

    unsigned int yUpper = pixel[1];
    if (yUpper != imageDimension[1] - 1)
      ++yUpper;

    unsigned int zLower = pixel[2];
    if (zLower != 0)
      --zLower;

    unsigned int zUpper = pixel[2];
    if (zUpper != imageDimension[2] - 1)
      ++zUpper;

    itk::Index<3> pixelXLower = pixel;
    pixelXLower[0] = xLower;
    itk::Index<3> pixelXUpper = pixel;
    pixelXUpper[0] = xUpper;

    itk::Index<3> pixelYLower = pixel;
    pixelYLower[1] = yLower;
    itk::Index<3> pixelYUpper = pixel;
    pixelYUpper[1] = yUpper;

    itk::Index<3> pixelZLower = pixel;
    pixelZLower[2] = zLower;
    itk::Index<3> pixelZUpper = pixel;
    pixelZUpper[2] = zUpper;

    if ((imagePixelWriter.GetPixelByIndex(pixelXLower) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if ((imagePixelWriter.GetPixelByIndex(pixelYLower) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if ((imagePixelWriter.GetPixelByIndex(pixelZLower) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if ((imagePixelWriter.GetPixelByIndex(pixelXUpper) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if ((imagePixelWriter.GetPixelByIndex(pixelYUpper) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if ((imagePixelWriter.GetPixelByIndex(pixelZUpper) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
  }
  return true;
}

void AblationUtils::CreateSafetyMarginInfluenceAreaOfPixel(itk::Index<3> &pixel,
                                                           mitk::Image::Pointer image,
                                                           double &margin,
                                                           mitk::Vector3D &imageDimension,
                                                           mitk::Vector3D &imageSpacing)
{
  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(margin / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(margin / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(margin / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           pixel,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (margin >= CalculateScalarDistance(pixel, actualIndex, imageSpacing))
          {
            if ((imagePixelWriter.GetPixelByIndex(actualIndex) & TUMOR_NOT_YET_ABLATED) != TUMOR_NOT_YET_ABLATED &&
                imagePixelWriter.GetPixelByIndex(actualIndex) != SAFETY_MARGIN)
            {
              imagePixelWriter.SetPixelByIndex(actualIndex, SAFETY_MARGIN);
            }
          }
        }
      }
    }
  }
}

double AblationUtils::CalculateRatioAblatedTissueOutsideTumorToAblatedTissueInsideTumor(itk::Index<3> &center,
                                                                                        mitk::Image::Pointer image,
                                                                                        double &radius,
                                                                                        mitk::Vector3D &imageDimension,
                                                                                        mitk::Vector3D &imageSpacing)
{
  double pixelsOutsideTumorAndSafetyMargin = 0.0;
  int totalNumberOfPixels = 0;

  if (image.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / imageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / imageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / imageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    CalculateUpperLowerXYZ(upperX,
                           lowerX,
                           upperY,
                           lowerY,
                           upperZ,
                           lowerZ,
                           pixelDirectionX,
                           pixelDirectionY,
                           pixelDirectionZ,
                           center,
                           imageDimension);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(center, actualIndex, imageSpacing))
          {
            if (imagePixelWriter.GetPixelByIndex(actualIndex) != TUMOR_NOT_YET_ABLATED &&
                imagePixelWriter.GetPixelByIndex(actualIndex) != SAFETY_MARGIN)
            {
              ++pixelsOutsideTumorAndSafetyMargin;
              ++totalNumberOfPixels;
            }
            else
            {
              ++totalNumberOfPixels;
            }
          }
        }
      }
    }
  }

  if (totalNumberOfPixels == 0)
  {
    // Prevent deviding by zero:
    totalNumberOfPixels = 1;
  }
  return pixelsOutsideTumorAndSafetyMargin / totalNumberOfPixels;
}

void AblationUtils::MoveCenterOfAblationZone(itk::Index<3> &center,
                                             mitk::Image::Pointer image,
                                             double &radius,
                                             mitk::Vector3D &imageDimension,
                                             mitk::Vector3D &imageSpacing)
{
  double distanceLowerX = 0.0;
  double distanceUpperX = 0.0;
  double distanceLowerY = 0.0;
  double distanceUpperY = 0.0;
  double distanceLowerZ = 0.0;
  double distanceUpperZ = 0.0;

  double distanceMoveUpX = 0.0;
  double distanceMoveDownX = 0.0;
  double distanceMoveUpY = 0.0;
  double distanceMoveDownY = 0.0;
  double distanceMoveUpZ = 0.0;
  double distanceMoveDownZ = 0.0;

  int resultMovingX = 0;
  int resultMovingY = 0;
  int resultMovingZ = 0;

  double radiusFactor = radius - 3;
  CalculateDistancesOfTumorBoundariesFromCenter(distanceLowerX,
                                                distanceUpperX,
                                                distanceLowerY,
                                                distanceUpperY,
                                                distanceLowerZ,
                                                distanceUpperZ,
                                                center,
                                                image,
                                                imageDimension,
                                                imageSpacing);

  MITK_WARN << "Distances lowerX " << distanceLowerX << " upperX " << distanceUpperX << " lowerY " << distanceLowerY
            << " upperY " << distanceUpperY << " lowerZ " << distanceLowerZ << " upperZ " << distanceUpperZ;
  if (distanceLowerX < radiusFactor)
  {
    distanceMoveUpX = radiusFactor - distanceLowerX;
    if ((distanceLowerX + distanceMoveUpX) > ((distanceLowerX + distanceUpperX) / 2))
    {
      distanceMoveUpX = ((distanceLowerX + distanceUpperX) / 2) - distanceLowerX;
    }
    distanceMoveUpX -= 2.0;
  }

  if (distanceUpperX < radiusFactor)
  {
    distanceMoveDownX = radiusFactor - distanceUpperX;
    if ((distanceUpperX + distanceMoveDownX) > ((distanceLowerX + distanceUpperX) / 2))
    {
      distanceMoveDownX = ((distanceLowerX + distanceUpperX) / 2) - distanceUpperX;
    }
    distanceMoveDownX -= 2.0;
  }

  if (distanceLowerY < radiusFactor)
  {
    distanceMoveUpY = radiusFactor - distanceLowerY;
    if ((distanceLowerY + distanceMoveUpY) > ((distanceLowerY + distanceUpperY) / 2))
    {
      distanceMoveUpY = ((distanceLowerY + distanceUpperY) / 2) - distanceLowerY;
    }
    distanceMoveUpY -= 2.0;
  }

  if (distanceUpperY < radiusFactor)
  {
    distanceMoveDownY = radiusFactor - distanceUpperY;
    if ((distanceUpperY + distanceMoveDownY) > ((distanceLowerY + distanceUpperY) / 2))
    {
      distanceMoveDownY = ((distanceLowerY + distanceUpperY) / 2) - distanceUpperY;
    }
    distanceMoveDownY -= 2.0;
  }

  if (distanceLowerZ < radiusFactor)
  {
    distanceMoveUpZ = radiusFactor - distanceLowerZ;
    if ((distanceLowerZ + distanceMoveUpZ) > ((distanceLowerZ + distanceUpperZ) / 2))
    {
      distanceMoveUpZ = ((distanceLowerZ + distanceUpperZ) / 2) - distanceLowerZ;
    }
    distanceMoveUpZ -= 2.0;
  }

  if (distanceUpperZ < radiusFactor)
  {
    distanceMoveDownZ = radiusFactor - distanceUpperZ;
    if ((distanceUpperZ + distanceMoveDownZ) > ((distanceLowerZ + distanceUpperZ) / 2))
    {
      distanceMoveDownZ = ((distanceLowerZ + distanceUpperZ) / 2) - distanceUpperZ;
    }
    distanceMoveDownZ -= 2.0;
  }

  resultMovingX = floor((distanceMoveUpX - distanceMoveDownX) / imageSpacing[0]);
  resultMovingY = floor((distanceMoveUpY - distanceMoveDownY) / imageSpacing[1]);
  resultMovingZ = floor((distanceMoveUpZ - distanceMoveDownZ) / imageSpacing[2]);

  MITK_INFO << "Moving ablation zone from: " << center[0] << "|" << center[1] << "|" << center[2];
  center[0] += resultMovingX;
  center[1] += resultMovingY;
  center[2] += resultMovingZ;
  MITK_INFO << "... to center: " << center[0] << "|" << center[1] << "|" << center[2];
}

int AblationUtils::CalculateTumorVolume(mitk::Image::Pointer image,
                                        mitk::Vector3D &imageSpacing,
                                        std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices)
{
  int volume = 0;
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    unsigned short pixelValue;
    for (int index = 0; index < tumorTissueSafetyMarginIndices.size(); ++index)
    {
      pixelValue = imagePixelWriter.GetPixelByIndex(tumorTissueSafetyMarginIndices.at(index));
      volume += (pixelValue & TUMOR_NOT_YET_ABLATED);
    }
  }
  double volumeFactor = imageSpacing[0] * imageSpacing[1] * imageSpacing[2];
  volume = (int)(volume * volumeFactor) / 1000;
  return volume;
}

int AblationUtils::CalculateSafetyMarginVolume(mitk::Image::Pointer image,
                                               mitk::Vector3D &imageSpacing,
                                               std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices)
{
  int volume = 0;
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    unsigned short pixelValue;
    for (int index = 0; index < tumorTissueSafetyMarginIndices.size(); ++index)
    {
      pixelValue = imagePixelWriter.GetPixelByIndex(tumorTissueSafetyMarginIndices.at(index));
      volume += ((pixelValue & SAFETY_MARGIN) / SAFETY_MARGIN);
    }
  }
  double volumeFactor = imageSpacing[0] * imageSpacing[1] * imageSpacing[2];
  volume = (int)(volume * volumeFactor) / 1000;
  return volume;
}

int AblationUtils::CalculateTotalAblationVolume(mitk::Image::Pointer image,
                                                mitk::Vector3D &imageSpacing,
                                                mitk::Vector3D &imageDimension)
{
  int volume = 0;
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
          pixelValue &= BIT_OPERATION_ELIMINATE_TUMOR_SAFETY_MARGIN;
          if (pixelValue != 0)
          {
            ++volume;
          }
        }
      }
    }
  }
  double volumeFactor = imageSpacing[0] * imageSpacing[1] * imageSpacing[2];
  volume = (int)(volume * volumeFactor) / 1000;
  return volume;
}

int AblationUtils::CalculateAblationVolumeAblatedMoreThanOneTime(mitk::Image::Pointer image,
                                                                 mitk::Vector3D &imageSpacing,
                                                                 mitk::Vector3D &imageDimension)
{
  int volume = 0;
  if (image.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
          pixelValue &= BIT_OPERATION_ELIMINATE_TUMOR_SAFETY_MARGIN;
          if (pixelValue != 0)
          {
            if ((pixelValue - ABLATION_VALUE) > 0)
            {
              ++volume;
            }
          }
        }
      }
    }
  }
  double volumeFactor = imageSpacing[0] * imageSpacing[1] * imageSpacing[2];
  volume = (int)(volume * volumeFactor) / 1000;
  return volume;
}
