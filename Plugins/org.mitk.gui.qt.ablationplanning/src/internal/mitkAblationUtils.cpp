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
#include "mitkAblationPlan.h"
#include "mitkProperties.h"
#include <cmath>
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

// rng

#if defined(_MSC_VER) // Visual studio
#define thread_local __declspec(thread)
#elif defined(__GCC__) // GCC
#define thread_local __thread
#endif

#include <random>
#include <thread>
#include <time.h>

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

int AblationUtils::intRand(const int &min, const int &max)
{
  int i;
  std::vector<int> v;

  static thread_local std::mt19937 *generator = nullptr;
  if (!generator)
    generator = new std::mt19937(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
  std::uniform_int_distribution<int> distribution(min, max);
  return distribution(*generator);
}

QString AblationUtils::FindAblationStartingPosition(mitk::Image::Pointer image,
                                                    std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices,
                                                    double ablationRadius,
                                                    double minRadius,
                                                    double maxRadius,
                                                    itk::Index<3> &tempAblationStartingPositionIndexCoordinates,
                                                    mitk::Point3D &tempAblationStartingPositionInWorldCoordinates,
                                                    double &tempAblationStartingRadius,
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

      double radius1 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        startingPosition1, image, imageSpacing, imageDimension, ablationRadius, minRadius, maxRadius);
      double radius2 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        startingPosition2, image, imageSpacing, imageDimension, ablationRadius, minRadius, maxRadius);
      double radius3 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        startingPosition3, image, imageSpacing, imageDimension, ablationRadius, minRadius, maxRadius);
      double radius4 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        startingPosition4, image, imageSpacing, imageDimension, ablationRadius, minRadius, maxRadius);
      double radius5 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        startingPosition5, image, imageSpacing, imageDimension, ablationRadius, minRadius, maxRadius);
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
      tempAblationStartingRadius = radiusVector.at(index);
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
                                            double radius,
                                            mitk::Vector3D &imageSpacing,
                                            mitk::Vector3D &imageDimension,
                                            std::vector<mitk::AblationZone> &tempAblationZones)
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
    tempAblationZones.push_back({center, radius});
  }
}

void AblationUtils::CalculateAblationVolume(itk::Index<3> &center,
                                            mitk::Image::Pointer image,
                                            double radius,
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
            // if (imagePixelWriter.GetPixelByIndex(actualIndex) == NO_TUMOR_ISSUE) other variant
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

double AblationUtils::GetPercentageOfVolumeInsideTumor(double &radius,
                                                       itk::Index<3> &centerOfVolume,
                                                       mitk::Image::Pointer image,
                                                       mitk::Vector3D &imageSpacing,
                                                       mitk::Vector3D &imageDimension)
{
  if (image.IsNull())
  {
    MITK_WARN << "No image present!";
    return 0.0;
  }

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
  int numberOfPixelsInsideTumor = 0;
  int numberOfPixelsOutsideTumor = 0;
  for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
  {
    for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
    {
      for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
      {
        if (radius >= CalculateScalarDistance(centerOfVolume, actualIndex, imageSpacing))
        {
          if (imagePixelWriter.GetPixelByIndex(actualIndex) == NO_TUMOR_ISSUE)
          {
            numberOfPixelsOutsideTumor++;
          }
          else if (imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED ||
                   imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN ||
                   imagePixelWriter.GetPixelByIndex(actualIndex) == ABLATION_VALUE)
          {
            numberOfPixelsInsideTumor++;
          }
        }
      }
    }
  }
  return (double)numberOfPixelsInsideTumor / (numberOfPixelsInsideTumor + numberOfPixelsOutsideTumor) * 100;
}

double AblationUtils::CalculateRadiusOfVolumeInsideTumorForGivenPoint(itk::Index<3> &point,
                                                                      mitk::Image::Pointer image,
                                                                      mitk::Vector3D &imageSpacing,
                                                                      mitk::Vector3D &imageDimension,
                                                                      const double startRadius,
                                                                      const double minRadius,
                                                                      const double maxRadius)
{
  {
    // double radius = startRadius;
    // while (CheckIfVolumeOfGivenRadiusIsTotallyInsideTumorTissueAndSafetyMargin(
    //         radius, point, image, imageSpacing, imageDimension) &&
    //       radius <= maxRadius)
    //{
    //  radius += 1;
    //}

    ///* Other variant: check percentage > 50 instead of whole zone inside
    // while (GetPercentageOfVolumeInsideTumor(radius, point, image, imageSpacing, imageDimension) > 50 && radius <=
    // maxRadius)
    //{
    //  radius += 1;
    //}
    //*/
    //// MITK_INFO << "Calculated max radius for given point " << point << " : " << radius;
    // return radius;

    double radius = startRadius;
    std::vector<double> radiusPercentage;
    bool maxRadiusReached{false};
    while (!maxRadiusReached)
    {
      radiusPercentage.push_back(
        GetPercentageOfTumorvolumeInsideZone(point, image, imageSpacing, imageDimension, radius));
      radius = radius + 1;
      if (radius >= maxRadius)
      {
        maxRadiusReached = true;
      }
    }

    for (int i = radiusPercentage.size() - 1; i >= 0; i--)
    {
      if (radiusPercentage.at(i) >= 0.7)
      {
        return startRadius + i;
      }
    }
    maxRadiusReached = false;
    std::vector<double> radiusNonAblatedPercentage;
    radius = startRadius;
    while (!maxRadiusReached)
    {
      radiusNonAblatedPercentage.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(point, image, imageSpacing, imageDimension, radius));
      radius += 1;
      if (radius >= maxRadius)
      {
        maxRadiusReached = true;
      }
    }
    int indexBestNonAblatedPercentage = 0;
    for (int k = 1; k < radiusNonAblatedPercentage.size(); k++)
    {
      if (radiusNonAblatedPercentage.at(indexBestNonAblatedPercentage) < radiusNonAblatedPercentage.at(k))
      {
        indexBestNonAblatedPercentage = k;
      }
    }
    radius = startRadius + indexBestNonAblatedPercentage;
    return radius;
  }
}

double AblationUtils::GetPercentageOfTumorvolumeInsideZone(itk::Index<3> &centerOfVolume,
                                                           mitk::Image::Pointer image,
                                                           mitk::Vector3D &imageSpacing,
                                                           mitk::Vector3D &imageDimension,
                                                           double radius)
{
  if (image.IsNull())
  {
    return 0;
  }

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

  mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
  itk::Index<3> actualIndex;
  unsigned short pixelValue;
  int zonePixels{0};
  int pixelsAblationNecessary{0};
  for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
  {
    for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
    {
      for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
      {
        if (radius >= CalculateScalarDistance(centerOfVolume, actualIndex, imageSpacing))
        {
          pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
          if (pixelValue >= TUMOR_NOT_YET_ABLATED)
          {
            pixelsAblationNecessary++;
          }
          zonePixels++;
        }
      }
    }
  }
  if (pixelsAblationNecessary == 0)
  {
    return 0;
  }
  return 1.0 / (double(zonePixels) / double(pixelsAblationNecessary));
}

double AblationUtils::GetPercentageOfNonAblatedTumorvolumeInsideZone(itk::Index<3> &centerOfVolume,
                                                                     mitk::Image::Pointer image,
                                                                     mitk::Vector3D &imageSpacing,
                                                                     mitk::Vector3D &imageDimension,
                                                                     double radius)
{
  if (image.IsNull())
  {
    return 0;
  }

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

  mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
  itk::Index<3> actualIndex;
  unsigned short pixelValue;
  int zonePixels{0};
  int pixelsAblationNecessary{0};
  for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
  {
    for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
    {
      for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
      {
        if (radius >= CalculateScalarDistance(centerOfVolume, actualIndex, imageSpacing))
        {
          pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
          if (pixelValue == TUMOR_NOT_YET_ABLATED || pixelValue == SAFETY_MARGIN)
          {
            pixelsAblationNecessary++;
          }
          zonePixels++;
        }
      }
    }
  }
  if (pixelsAblationNecessary == 0)
  {
    return 0;
  }
  return double(pixelsAblationNecessary) / double(zonePixels);
}

int AblationUtils::GetNumberOfAblatedPoints(itk::Index<3> &centerOfVolume,
                                            mitk::Image::Pointer image,
                                            mitk::Vector3D &imageSpacing,
                                            mitk::Vector3D &imageDimension,
                                            double radius)
{
  if (image.IsNull())
  {
    return 0;
  }

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

  mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
  itk::Index<3> actualIndex;
  unsigned short pixelValue;
  int zonePixels{0};
  int pixelAblatedOnce{0};
  for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
  {
    for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
    {
      for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
      {
        if (radius >= CalculateScalarDistance(centerOfVolume, actualIndex, imageSpacing))
        {
          pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
          if (pixelValue >= 255)
          {
            pixelValue -= 255;
          }
          if (pixelValue < 4 || pixelValue > 1)
          {
            pixelAblatedOnce++;
          }
        }
      }
    }
  }
  return pixelAblatedOnce;
}

double AblationUtils::CheckImageForNonAblatedTissueInPercentage(mitk::Image::Pointer image,
                                                                mitk::Vector3D &imageDimension)
{
  int numberOfNonAblatedPixels = 0;
  int numberOfTumorPixels = 0;
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
            numberOfTumorPixels++;
            numberOfNonAblatedPixels++;
          }
          else if (imagePixelWriter.GetPixelByIndex(actualIndex) == ABLATION_VALUE)
          {
            numberOfTumorPixels++;
          }
        }
      }
    }
  }
  return ((double)numberOfNonAblatedPixels / numberOfTumorPixels) * 100;
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

bool AblationUtils::CheckForNonAblatedTumorTissueWithSafetyMargin(std::vector<itk::Index<3>> &indices,
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
          if (imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN &&
              imagePixelWriter.GetPixelByIndex(actualIndex) == ABLATION_VALUE)
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

void AblationUtils::DetectNotNeededAblationVolume(mitk::AblationPlan::Pointer plan,
                                                  mitk::Image::Pointer image,
                                                  mitk::Vector3D &imageDimension,
                                                  mitk::Vector3D &imageSpacing)
{
  std::vector<int> indicesRemoved;
  for (int index = 0; index < plan->GetNumberOfZones(); ++index)
  {
    if (!CheckIfAblationVolumeIsNeeded(plan->GetAblationZone(index)->indexCenter,
                                       image,
                                       plan->GetAblationZone(index)->radius,
                                       imageDimension,
                                       imageSpacing))
    {
      RemoveAblationVolume(plan->GetAblationZone(index)->indexCenter,
                           image,
                           plan->GetAblationZone(index)->radius,
                           imageDimension,
                           imageSpacing);
      indicesRemoved.push_back(index);
    }
  }
  for (int index = indicesRemoved.size() - 1; index >= 0; --index)
  {
    MITK_INFO << "Radius ablation Zone not needed: " << plan->GetAblationZone(indicesRemoved.at(index))->radius;
    plan->RemoveAblationZone(indicesRemoved.at(index));
    MITK_INFO << "Removed Ablation zone at index position: " << indicesRemoved.at(index);
  }
}

void AblationUtils::RemoveNotNeededAblationZones(mitk::AblationPlan::Pointer plan,
                                                 mitk::Image::Pointer image,
                                                 mitk::Vector3D &imageDimension,
                                                 mitk::Vector3D &imageSpacing,
                                                 std::vector<itk::Index<3>> &m_TumorTissueSafetyMarginIndices,
                                                 double m_ToleranceNonAblatedTumorSafetyMarginVolume)
{
  if (image.IsNotNull())
  {
    double volumeTumor = AblationUtils::CalculateTumorVolume(image, imageSpacing, m_TumorTissueSafetyMarginIndices);
    while (true)
    {
      std::vector<double> overlappingVolume;
      double factorNonAblatedVolume = AblationUtils::CheckImageForNonAblatedTissueInPercentage(
        plan->GetSegmentationImage(), plan->GetImageDimension());
      {
        mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
        for (int i = 0; i < plan->GetNumberOfZones(); i++)
        {
          overlappingVolume.push_back(0);

          unsigned int pixelDirectionX = floor(plan->GetAblationZone(i)->radius / imageSpacing[0]);
          unsigned int pixelDirectionY = floor(plan->GetAblationZone(i)->radius / imageSpacing[1]);
          unsigned int pixelDirectionZ = floor(plan->GetAblationZone(i)->radius / imageSpacing[2]);

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
                                 plan->GetAblationZone(i)->indexCenter,
                                 imageDimension);
          for (int index0 = lowerX; index0 <= upperX; index0++)
          {
            for (int index1 = lowerY; index1 <= upperY; index1++)
            {
              for (int index2 = lowerZ; index2 <= upperZ; index2++)
              {
                itk::Index<3> actualIndex = {index0, index1, index2};
                if (AblationUtils::CalculateScalarDistance(plan->GetAblationZone(i)->indexCenter,
                                                           actualIndex,
                                                           imageSpacing) <= plan->GetAblationZone(i)->radius)
                {
                  int pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
                  if (pixelValue >= SAFETY_MARGIN)
                  {
                    pixelValue -= 255;
                  }
                  if (pixelValue - ABLATION_VALUE <= ABLATION_VALUE && pixelValue % 2 == TUMOR_NOT_YET_ABLATED)
                  {
                    overlappingVolume[i]++;
                  }
                }
              }
            }
          }
        }
      }
      int zoneToDelete{0};
      for (int i = 1; i < overlappingVolume.size(); i++) // find smallest
      {
        if (overlappingVolume[i] < overlappingVolume[zoneToDelete])
        {
          zoneToDelete = i;
        }
      }
      double imageSpacingFactor = imageSpacing[0] * imageSpacing[1] * imageSpacing[2];
      overlappingVolume[zoneToDelete] = overlappingVolume[zoneToDelete] * imageSpacingFactor / 1000.0;
      double overlappingFactor = (overlappingVolume[zoneToDelete] / volumeTumor) * 100;
      if (overlappingFactor + factorNonAblatedVolume < m_ToleranceNonAblatedTumorSafetyMarginVolume * 100)
      {
        RemoveAblationVolume(plan->GetAblationZone(zoneToDelete)->indexCenter,
                             image,
                             plan->GetAblationZone(zoneToDelete)->radius,
                             imageDimension,
                             imageSpacing);
        plan->RemoveAblationZone(zoneToDelete);

        factorNonAblatedVolume += overlappingFactor;
      }
      else
      {
        return;
      }
    }
    return;
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

    while (currentRadius > minRadius)
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
                return currentRadius;
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
    return;
  }
  return;
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

mitk::AblationZone AblationUtils::SearchNextAblationCenter(std::vector<itk::Index<3>> &tumorSafetyMarginPixels,
                                                           std::vector<itk::Index<3>> &unchangedTumorSafetyMarginPixels,
                                                           mitk::Image::Pointer image,
                                                           double &radius,
                                                           double &minRadius,
                                                           double &maxRadius,
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
      int randomIndex6 = rand() % unchangedTumorSafetyMarginPixels.size();
      int randomIndex7 = rand() % unchangedTumorSafetyMarginPixels.size();
      int randomIndex8 = rand() % unchangedTumorSafetyMarginPixels.size();
      int randomIndex9 = rand() % unchangedTumorSafetyMarginPixels.size();
      int randomIndex10 = rand() % unchangedTumorSafetyMarginPixels.size();

      itk::Index<3> position1 = tumorSafetyMarginPixels.at(randomIndex1);
      itk::Index<3> position2 = tumorSafetyMarginPixels.at(randomIndex2);
      itk::Index<3> position3 = tumorSafetyMarginPixels.at(randomIndex3);
      itk::Index<3> position4 = tumorSafetyMarginPixels.at(randomIndex4);
      itk::Index<3> position5 = unchangedTumorSafetyMarginPixels.at(randomIndex5);
      itk::Index<3> position6 = unchangedTumorSafetyMarginPixels.at(randomIndex6);
      itk::Index<3> position7 = unchangedTumorSafetyMarginPixels.at(randomIndex7);
      itk::Index<3> position8 = unchangedTumorSafetyMarginPixels.at(randomIndex8);
      itk::Index<3> position9 = unchangedTumorSafetyMarginPixels.at(randomIndex9);
      itk::Index<3> position10 = unchangedTumorSafetyMarginPixels.at(randomIndex10);

      std::vector<itk::Index<3>> positions;
      positions.push_back(position1);
      positions.push_back(position2);
      positions.push_back(position3);
      positions.push_back(position4);
      positions.push_back(position5);
      positions.push_back(position6);
      positions.push_back(position7);
      positions.push_back(position8);
      positions.push_back(position9);
      positions.push_back(position10);

      double radius1 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position1, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius2 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position2, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius3 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position3, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius4 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position4, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius5 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position5, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius6 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position6, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius7 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position7, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius8 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position8, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius9 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position9, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      double radius10 = CalculateRadiusOfVolumeInsideTumorForGivenPoint(
        position10, image, imageSpacing, imageDimension, radius, minRadius, maxRadius);
      std::vector<double> radiusVector;
      std::vector<int>::iterator result;
      std::vector<double> percentageNonAblatedTumorVector;
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position1, image, imageSpacing, imageDimension, radius1));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position2, image, imageSpacing, imageDimension, radius2));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position3, image, imageSpacing, imageDimension, radius3));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position4, image, imageSpacing, imageDimension, radius4));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position5, image, imageSpacing, imageDimension, radius5));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position6, image, imageSpacing, imageDimension, radius6));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position7, image, imageSpacing, imageDimension, radius7));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position8, image, imageSpacing, imageDimension, radius8));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position9, image, imageSpacing, imageDimension, radius9));
      percentageNonAblatedTumorVector.push_back(
        GetPercentageOfNonAblatedTumorvolumeInsideZone(position10, image, imageSpacing, imageDimension, radius10));
          //GetPercentageOfNonAblatedTumorvolumeInsideZone
      std::vector<double> percentageTumorVector;
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position1, image, imageSpacing, imageDimension, radius1));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position2, image, imageSpacing, imageDimension, radius2));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position3, image, imageSpacing, imageDimension, radius3));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position4, image, imageSpacing, imageDimension, radius4));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position5, image, imageSpacing, imageDimension, radius5));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position6, image, imageSpacing, imageDimension, radius6));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position7, image, imageSpacing, imageDimension, radius7));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position8, image, imageSpacing, imageDimension, radius8));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position9, image, imageSpacing, imageDimension, radius9));
      percentageTumorVector.push_back(
        GetPercentageOfTumorvolumeInsideZone(position10, image, imageSpacing, imageDimension, radius10));

      radiusVector.push_back(radius1);
      radiusVector.push_back(radius2);
      radiusVector.push_back(radius3);
      radiusVector.push_back(radius4);
      radiusVector.push_back(radius5);
      radiusVector.push_back(radius6);
      radiusVector.push_back(radius7);
      radiusVector.push_back(radius8);
      radiusVector.push_back(radius9);
      radiusVector.push_back(radius10);

      // std::vector<int> ablatedPoints;

      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position1, image, imageSpacing, imageDimension, radius1));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position2, image, imageSpacing, imageDimension, radius2));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position3, image, imageSpacing, imageDimension, radius3));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position4, image, imageSpacing, imageDimension, radius4));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position5, image, imageSpacing, imageDimension, radius5));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position6, image, imageSpacing, imageDimension, radius6));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position7, image, imageSpacing, imageDimension, radius7));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position8, image, imageSpacing, imageDimension, radius8));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position9, image, imageSpacing, imageDimension, radius9));
      // ablatedPoints.push_back(GetNumberOfAblatedPoints(position10, image, imageSpacing, imageDimension, radius10));

      int index = {0};
      for (int i = 1; i < radiusVector.size(); i++)
      {
        if (percentageNonAblatedTumorVector.at(i) > percentageNonAblatedTumorVector.at(index) && percentageTumorVector.at(i) > 0.5)
        {
          index = i;
        }
      }
      if (percentageNonAblatedTumorVector.at(index) > 0.2 && percentageTumorVector.at(index) >= 0.5)
      {
        return {positions.at(index), radiusVector.at(index)};
      }
      else
      {
        ++iteration;
      }
      if (iteration == 19)
      {
        return {positions.at(index), radiusVector.at(index)};
      }
    }
  }
  return {itk::Index<3>(), 0.0};
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

void AblationUtils::MoveCenterTowardsCenter(itk::Index<3> &center,
                                            mitk::Image::Pointer image,
                                            mitk::Vector3D &imageDimension,
                                            mitk::Vector3D &imageSpacing,
                                            double startRadius)
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

  double radiusFactor = startRadius;
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

  // MITK_WARN << "Distances lowerX " << distanceLowerX << " upperX " << distanceUpperX << " lowerY " <<
  // distanceLowerY
  //          << " upperY " << distanceUpperY << " lowerZ " << distanceLowerZ << " upperZ " << distanceUpperZ;
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

  // MITK_INFO << "Moving ablation zone from: " << center[0] << "|" << center[1] << "|" << center[2];
  center[0] += resultMovingX;
  center[1] += resultMovingY;
  center[2] += resultMovingZ;
  // MITK_INFO << "... to center: " << center[0] << "|" << center[1] << "|" << center[2];
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

  // MITK_WARN << "Distances lowerX " << distanceLowerX << " upperX " << distanceUpperX << " lowerY " <<
  // distanceLowerY
  //          << " upperY " << distanceUpperY << " lowerZ " << distanceLowerZ << " upperZ " << distanceUpperZ;
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

  // MITK_INFO << "Moving ablation zone from: " << center[0] << "|" << center[1] << "|" << center[2];
  center[0] += resultMovingX;
  center[1] += resultMovingY;
  center[2] += resultMovingZ;
  // MITK_INFO << "... to center: " << center[0] << "|" << center[1] << "|" << center[2];
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
    mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          unsigned short pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
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
  int volume{0};
  if (image.IsNotNull())
  {
    mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0] += 1)
        {
          unsigned short pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
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

void AblationUtils::ComputeStatistics(mitk::AblationPlan::Pointer plan,
                                      std::vector<itk::Index<3>> tumorTissueSafetyMarginIndices,
                                      double factorMaxNonAblatedVolume)
{
  mitk::AblationPlan::AblationPlanStatistics s;
  s.tumorVolume = AblationUtils::CalculateTumorVolume(
    plan->GetSegmentationImage(), plan->GetImageSpacing(), tumorTissueSafetyMarginIndices);
  s.safetyMarginVolume = AblationUtils::CalculateSafetyMarginVolume(
    plan->GetSegmentationImage(), plan->GetImageSpacing(), tumorTissueSafetyMarginIndices);
  s.tumorAndSafetyMarginVolume = s.tumorVolume + s.safetyMarginVolume;
  s.totalAblationVolume = AblationUtils::CalculateTotalAblationVolume(
    plan->GetSegmentationImage(), plan->GetImageSpacing(), plan->GetImageDimension());
  s.ablationVolumeAblatedMoreThanOneTime = AblationUtils::CalculateAblationVolumeAblatedMoreThanOneTime(
    plan->GetSegmentationImage(), plan->GetImageSpacing(), plan->GetImageDimension());

  s.factorOverlappingAblationZones = ((double)s.ablationVolumeAblatedMoreThanOneTime / s.totalAblationVolume) * 100;

  s.factorAblatedVolumeOutsideSafetyMargin =
    ((s.totalAblationVolume - s.tumorAndSafetyMarginVolume) / (double)s.totalAblationVolume) * 100;
  s.factorNonAblatedVolume =
    AblationUtils::CheckImageForNonAblatedTissueInPercentage(plan->GetSegmentationImage(), plan->GetImageDimension());
  s.factorMaxNonAblatedVolume = factorMaxNonAblatedVolume;
  //AblationUtils::FindAgglomerations(
  //  plan->GetSegmentationImage(), plan->GetImageSpacing(), plan->GetImageDimension());
  plan->SetStatistics(s);
}

bool AblationUtils::CheckIfVolumeMostlyInsideTumorAndSafetymarginTissue(
  itk::Index<3> &zoneCenter,
  std::vector<itk::Index<3>> &tumorSafetyMarginPixels,
  mitk::Image::Pointer &image,
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
                           zoneCenter,
                           imageDimension);

    mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
    itk::Index<3> actualIndex;
    int safetyAndTumorPixels{0};
    int totalPixels{0};
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (radius >= CalculateScalarDistance(zoneCenter, actualIndex, imageSpacing))
          {
            totalPixels++;
            int pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
            if (pixelValue > 255)
            {
              pixelValue -= 255;
            }
            if (pixelValue % 2 != 1)
            {
              safetyAndTumorPixels++;
            }
          }
        }
      }
    }
    if (double(safetyAndTumorPixels) / double(totalPixels) > 0.0)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

std::vector<std::vector<itk::Index<3>>> AblationUtils::FindAgglomerations(mitk::Image::Pointer image,
                                                                          mitk::Vector3D &imageSpacing,
                                                                          mitk::Vector3D &imageDimension)
{
  MITK_INFO << "Entered FindAgglomerations()";
  itk::Index<3> actualIndex;
  unsigned short pixelValue;
  std::vector<std::vector<itk::Index<3>>> agglomerationList;
  if (image.IsNotNull())
  {
    for (actualIndex[2] = 0; actualIndex[2] < imageDimension[2]; actualIndex[2]++)
    {
      for (actualIndex[1] = 0; actualIndex[1] < imageDimension[1]; actualIndex[1]++)
      {
        for (actualIndex[0] = 0; actualIndex[0] < imageDimension[0]; actualIndex[0]++)
        {
          {
            mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
            pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
          }
          if (pixelValue == TUMOR_NOT_YET_ABLATED || pixelValue == SAFETY_MARGIN)
          {
            if (!CheckIfPixelIsElementOfAgglomerationList(agglomerationList, actualIndex))
            {
              std::vector<itk::Index<3>> tempAgglomeration;
              tempAgglomeration = FindFullAgglomeration(image, actualIndex);
              agglomerationList.push_back(tempAgglomeration);
            }
          }
        } // pixelFinder
      }
    }
  }
  MITK_INFO << "Number of Agglomerations found:" << agglomerationList.size();
  for (size_t i = 0; i < agglomerationList.size(); i++)
  {
    MITK_INFO << "Liste Nr.: " << i << ", " << agglomerationList.at(i).size();
  }
  return agglomerationList;
}


std::vector<itk::Index<3>> AblationUtils::FindFullAgglomeration(mitk::Image::Pointer image,
                                                                itk::Index<3> &startingIndex)
{
  std::vector<itk::Index<3>> tempAgglomeration;
  tempAgglomeration.push_back(startingIndex);
  for (int i = 0; i < tempAgglomeration.size(); i++)
  {
    startingIndex = tempAgglomeration.at(i);
    itk::Index<3> actualIndex;
    unsigned short pixelValue;
    mitk::ImagePixelReadAccessor<unsigned short, 3> imagePixelReader(image);
    actualIndex = startingIndex;
    for (actualIndex[2] -= 1; actualIndex[2] <= startingIndex[2] + 1; actualIndex[2]++)
    {
      for (actualIndex[1] -= 1; actualIndex[1] <= startingIndex[1] + 1; actualIndex[1]++)
      {
        for (actualIndex[0] -= 1; actualIndex[0] <= startingIndex[0] + 1; actualIndex[0]++)
        {
          pixelValue = imagePixelReader.GetPixelByIndex(actualIndex);
          if (pixelValue == TUMOR_NOT_YET_ABLATED || pixelValue == SAFETY_MARGIN)
          {
            if (!CheckIfPixelIsElementOfAgglomeration(tempAgglomeration, actualIndex))
            {
              tempAgglomeration.push_back(actualIndex);
            }
          }
        }
      }
    }
  }
  return tempAgglomeration;
}

bool AblationUtils::CheckIfPixelIsElementOfAgglomeration(std::vector<itk::Index<3>> &agglomeration, itk::Index<3> index)
{
  for (size_t i = 0; i < agglomeration.size(); i++)
  {
    if (agglomeration.at(i) == index)
    {
      return true;
    }
  }
  return false;
}

bool AblationUtils::CheckIfPixelIsElementOfAgglomerationList(std::vector<std::vector<itk::Index<3>>> &agglomerationList,
                                                             itk::Index<3> pixel)
{
  for (int i = 0; i < agglomerationList.size(); i++)
  {
    for (int j = 0; j < agglomerationList.at(i).size(); j++)
    {
      if (agglomerationList.at(i).at(j) == pixel)
      {
        return true;
      }
    }
  }
  return false;
}

void AblationUtils::SetSolutionValueStatistics(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans)
{
  SetMinMaxAblationZoneNumber(AllFoundPlans);
  SetMinMaxOverlapVolume(AllFoundPlans);
  SetMinMaxVolumeOutsideFactor(AllFoundPlans);
  // SetMinMaxAgglomerations(AllFoundPlans); //still need to think of a concept
  return;
}

void AblationUtils::SetMinMaxAblationZoneNumber(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans)
{
  for (int i = 0; i < AllFoundPlans.size(); i++)
  {
    if (AllFoundPlans.at(i).GetPointer()->GetNumberOfZones() >
        AllFoundPlans.at(i).GetPointer()->GetMaxAblationZoneNumber())
    { // set max number of Zones
      for (int j = 0; j < AllFoundPlans.size(); j++)
      {
        AllFoundPlans.at(j).GetPointer()->SetMaxAblationZoneNumber(
          AllFoundPlans.at(i).GetPointer()->GetNumberOfZones());
      }
    } // set min number of Zones
    if (AllFoundPlans.at(i).GetPointer()->GetNumberOfZones() <
        AllFoundPlans.at(i).GetPointer()->GetMinAblationZoneNumber())
    {
      for (int j = 0; j < AllFoundPlans.size(); j++)
      {
        AllFoundPlans.at(j).GetPointer()->SetMinAblationZoneNumber(
          AllFoundPlans.at(i).GetPointer()->GetNumberOfZones());
      }
    }
  }
  return;
}

void AblationUtils::SetMinMaxOverlapVolume(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans)
{
  for (int i = 0; i < AllFoundPlans.size(); i++)
  {
    AllFoundPlans.at(i).GetPointer()->SetMinOverlappingZonesFactor(
      AllFoundPlans.at(0)->GetStatistics().factorOverlappingAblationZones);
  }
  for (int i = 0; i < AllFoundPlans.size(); i++)
  {
    if (AllFoundPlans.at(i)->GetStatistics().factorOverlappingAblationZones >
        AllFoundPlans.at(i).GetPointer()->GetMaxOverlappingZonesFactor())
    {
      // set max factor of overlapping Zones
      for (int j = 0; j < AllFoundPlans.size(); j++)
      {
        AllFoundPlans.at(j).GetPointer()->SetMaxOverlappingZonesFactor(
          AllFoundPlans.at(i)->GetStatistics().factorOverlappingAblationZones);
      }
    } // set min number of overlapping Zones
    if (AllFoundPlans.at(i)->GetStatistics().factorOverlappingAblationZones <
        AllFoundPlans.at(i).GetPointer()->GetMinOverlappingZonesFactor())
    {
      for (int j = 0; j < AllFoundPlans.size(); j++)
      {
        AllFoundPlans.at(j).GetPointer()->SetMinOverlappingZonesFactor(
          AllFoundPlans.at(i)->GetStatistics().factorOverlappingAblationZones);
      }
    }
  }
  return;
}

void AblationUtils::SetMinMaxVolumeOutsideFactor(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans)
{
  for (int i = 0; i < AllFoundPlans.size(); i++)
  {
    if (AllFoundPlans.at(i)->GetStatistics().factorAblatedVolumeOutsideSafetyMargin >
        AllFoundPlans.at(i).GetPointer()->GetMaxVolumeOutsideFactor())
    {
      // set max factor of volume outside
      for (int j = 0; j < AllFoundPlans.size(); j++)
      {
        AllFoundPlans.at(j).GetPointer()->SetMaxVolumeOutsideFactor(
          AllFoundPlans.at(i)->GetStatistics().factorAblatedVolumeOutsideSafetyMargin);
      }
    } // set min number of overlapping Zones
    if (AllFoundPlans.at(i)->GetStatistics().factorAblatedVolumeOutsideSafetyMargin <
        AllFoundPlans.at(i).GetPointer()->GetMinVolumeOutsideFactor())
    {
      for (int j = 0; j < AllFoundPlans.size(); j++)
      {
        AllFoundPlans.at(j).GetPointer()->SetMinVolumeOutsideFactor(
          AllFoundPlans.at(i)->GetStatistics().factorAblatedVolumeOutsideSafetyMargin);
      }
    }
  }
  return;
}