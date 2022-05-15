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

#ifndef ABLATIONUTILS_H
#define ABLATIONUTILS_H

#include "mitkAblationPlan.h"
#include "mitkAblationZone.h"
#include <QmitkAbstractView.h>
#include <berryISelectionListener.h>
#include <mitkImage.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateOr.h>

/**
  \brief AblationUtils static class for doing ablation calculations

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \ingroup ${plugin_target}_internal
*/

class AblationUtils
{
public:
  static void FillVectorContainingIndicesOfTumorTissueSafetyMargin(
    mitk::Image::Pointer image,
    mitk::Vector3D &imageDimension,
    std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices);

  static void ComputeStatistics(mitk::AblationPlan::Pointer plan,
                                std::vector<itk::Index<3>> tumorTissueSafetyMarginIndices,
                                double factorMaxNonAblatedTumor);

  static std::vector<itk::Index<3>> FillVectorContainingIndicesOfTumorTissueOnly(mitk::Image::Pointer image,
                                                                                 mitk::Vector3D &imageDimension);

  /** Searches for a new starting point for the first ablation sphere.
   *  @param[in] image  image on which the ablation is planned (needed for geometry calculations)
   *  @param[in] tumorTissueSafetyMarginIndices   Indices of tumor (+ safety margin) that are not ablated yet
   *  @param[in] ablationRadius   Desired radius of the ablation
   *  @param[out] tempAblationStartingPositionIndexCoordinates    Returns the corrdinates of the best starting point
   * that was found
   *  @param[out] tempAblationStartingPositionInWorldCoordinates  Returns the corrdinates of the best starting point
   * that was found
   *  @param[out] tempAblationStartingRadius                      Returns the radius of the starting ablation volume
   *  @param[in]  imageDimension
   *  @param[in]  imageSpacing
   */
  static QString FindAblationStartingPosition(mitk::Image::Pointer image,
                                              std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices,
                                              double ablationRadius,
                                              double minRadius,
                                              double maxRadius,
                                              itk::Index<3> &tempAblationStartingPositionIndexCoordinates,
                                              mitk::Point3D &tempAblationStartingPositionInWorldCoordinates,
                                              double &tempAblationStartingRadius,
                                              mitk::Vector3D &imageDimension,
                                              mitk::Vector3D &imageSpacing);

  static double CalculateScalarDistance(itk::Index<3> &point1, itk::Index<3> &point2, mitk::Vector3D &imageSpacing);

  static void CalculateAblationVolume(itk::Index<3> &center,
                                      mitk::Image::Pointer image,
                                      double radius,
                                      mitk::Vector3D &imageSpacing,
                                      mitk::Vector3D &imageDimension,
                                      std::vector<mitk::AblationZone> &tempAblationZones);

  static void CalculateAblationVolume(itk::Index<3> &center,
                                      mitk::Image::Pointer image,
                                      double radius,
                                      mitk::Vector3D &imageSpacing,
                                      mitk::Vector3D &imageDimension);

  static bool CheckVolumeForNonAblatedTissue(itk::Index<3> &centerOfVolume,
                                             mitk::Image::Pointer image,
                                             double &radius,
                                             mitk::Vector3D &imageSpacing,
                                             mitk::Vector3D &imageDimension);

  /** @return Returns the percentage of a sphere with the given center and radius that is inside the tumor volume of the
   * given image */
  static double GetPercentageOfVolumeInsideTumor(double &radius,
                                                 itk::Index<3> &centerOfVolume,
                                                 mitk::Image::Pointer image,
                                                 mitk::Vector3D &imageSpacing,
                                                 mitk::Vector3D &imageDimension);

  /** Checks if a sphere with the given center and radius is totally inside the tumor volume of the given image */
  static bool CheckIfVolumeOfGivenRadiusIsTotallyInsideTumorTissueAndSafetyMargin(double &radius,
                                                                                  itk::Index<3> &centerOfVolume,
                                                                                  mitk::Image::Pointer image,
                                                                                  mitk::Vector3D &imageSpacing,
                                                                                  mitk::Vector3D &imageDimension);

  static double CalculateRadiusOfVolumeInsideTumorForGivenPoint(itk::Index<3> &point,
                                                                mitk::Image::Pointer image,
                                                                mitk::Vector3D &imageSpacing,
                                                                mitk::Vector3D &imageDimension,
                                                                double startRadius,
                                                                double minRadius,
                                                                double maxRadius);

  /** @returns Returns the percentage of non ablated tumor tissue (based on voxels) */
  static double CheckImageForNonAblatedTissueInPercentage(mitk::Image::Pointer image, mitk::Vector3D &imageDimension);

  static bool CheckImageForNonAblatedTissue(mitk::Image::Pointer image, mitk::Vector3D &imageDimension);

  static bool CheckForNonAblatedTumorTissueWithoutSafetyMargin(std::vector<itk::Index<3>> &indices,
                                                               mitk::Image::Pointer image,
                                                               mitk::Vector3D &imageDimension);

  static bool CheckForNonAblatedTumorTissueWithSafetyMargin(std::vector<itk::Index<3>> &indices,
                                                            mitk::Image::Pointer image,
                                                            mitk::Vector3D &imageDimension);

  static void CalculateUpperLowerXYZ(unsigned int &upperX,
                                     unsigned int &lowerX,
                                     unsigned int &upperY,
                                     unsigned int &lowerY,
                                     unsigned int &upperZ,
                                     unsigned int &lowerZ,
                                     unsigned int &pixelDirectionX,
                                     unsigned int &pixelDirectionY,
                                     unsigned int &pixelDirectionZ,
                                     itk::Index<3> &center,
                                     mitk::Vector3D &imageDimension);

  static void CalculateDistancesOfTumorBoundariesFromCenter(double &distanceLowerX,
                                                            double &distanceUpperX,
                                                            double &distanceLowerY,
                                                            double &distanceUpperY,
                                                            double &distanceLowerZ,
                                                            double &distanceUpperZ,
                                                            itk::Index<3> &center,
                                                            mitk::Image::Pointer image,
                                                            mitk::Vector3D &imageDimension,
                                                            mitk::Vector3D &imageSpacing);

  static void DetectNotNeededAblationVolume(mitk::AblationPlan::Pointer plan,
                                            mitk::Image::Pointer image,
                                            mitk::Vector3D &imageDimension,
                                            mitk::Vector3D &imageSpacing);

  static void RemoveNotNeededAblationZones(mitk::AblationPlan::Pointer plan,
                                           mitk::Image::Pointer image,
                                           mitk::Vector3D &imageDimension,
                                           mitk::Vector3D &imageSpacing,
                                           std::vector<itk::Index<3>> &m_TumorTissueSafetyMarginIndices,
                                           double m_ToleranceNonAblatedTumorSafetyMarginVolume);

  static double FindMinimalAblationRadius(itk::Index<3> &center,
                                          mitk::Image::Pointer image,
                                          double &maxRadius,
                                          double &minRadius,
                                          mitk::Vector3D &imageDimension,
                                          mitk::Vector3D &imageSpacing);

  static bool CheckIfAblationVolumeIsNeeded(itk::Index<3> &center,
                                            mitk::Image::Pointer image,
                                            double &radius,
                                            mitk::Vector3D &imageDimension,
                                            mitk::Vector3D &imageSpacing);

  static void RemoveAblationVolume(itk::Index<3> &center,
                                   mitk::Image::Pointer image,
                                   double &radius,
                                   mitk::Vector3D &imageDimension,
                                   mitk::Vector3D &imageSpacing);

  static void RemoveAblatedPixelsFromGivenVector(itk::Index<3> &center,
                                                 std::vector<itk::Index<3>> &tumorSafetyMarginPixels,
                                                 mitk::Image::Pointer image,
                                                 double &radius,
                                                 mitk::Vector3D &imageDimension,
                                                 mitk::Vector3D &imageSpacing);

  static mitk::AblationZone SearchNextAblationCenter(std::vector<itk::Index<3>> &tumorSafetyMarginPixels,
                                                     std::vector<itk::Index<3>> &unchangedTumorSafetyMarginPixels,
                                                     mitk::Image::Pointer image,
                                                     double &radius,
                                                     double &minRadius,
                                                     double &maxRadius,
                                                     mitk::Vector3D &imageDimension,
                                                     mitk::Vector3D &imageSpacing);

  static void ResetSegmentationImage(mitk::Image::Pointer image, mitk::Vector3D &imageDimension);

  static void ResetSafetyMargin(mitk::Image::Pointer image, mitk::Vector3D &imageDimension);

  static bool CheckAllVonNeumannNeighbourPixelsAreTumorTissue(itk::Index<3> &pixel,
                                                              mitk::Image::Pointer image,
                                                              mitk::Vector3D &imageDimension);

  static void CreateSafetyMarginInfluenceAreaOfPixel(itk::Index<3> &pixel,
                                                     mitk::Image::Pointer image,
                                                     double &margin,
                                                     mitk::Vector3D &imageDimension,
                                                     mitk::Vector3D &imageSpacing);

  static double CalculateRatioAblatedTissueOutsideTumorToAblatedTissueInsideTumor(itk::Index<3> &center,
                                                                                  mitk::Image::Pointer image,
                                                                                  double &radius,
                                                                                  mitk::Vector3D &imageDimension,
                                                                                  mitk::Vector3D &imageSpacing);

  static void MoveCenterTowardsCenter(itk::Index<3> &center,
                                      mitk::Image::Pointer image,
                                      mitk::Vector3D &imageDimension,
                                      mitk::Vector3D &imageSpacing,
                                      double startRadius);

  static void MoveCenterOfAblationZone(itk::Index<3> &center,
                                       mitk::Image::Pointer image,
                                       double &radius,
                                       mitk::Vector3D &imageDimension,
                                       mitk::Vector3D &imageSpacing);

  static int CalculateTumorVolume(mitk::Image::Pointer image,
                                  mitk::Vector3D &imageSpacing,
                                  std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices);

  static int CalculateSafetyMarginVolume(mitk::Image::Pointer image,
                                         mitk::Vector3D &imageSpacing,
                                         std::vector<itk::Index<3>> &tumorTissueSafetyMarginIndices);

  static int CalculateTotalAblationVolume(mitk::Image::Pointer image,
                                          mitk::Vector3D &imageSpacing,
                                          mitk::Vector3D &imageDimension);

  static int CalculateAblationVolumeAblatedMoreThanOneTime(mitk::Image::Pointer image,
                                                           mitk::Vector3D &imageSpacing,
                                                           mitk::Vector3D &imageDimension);

  static bool CheckIfVolumeMostlyInsideTumorAndSafetymarginTissue(itk::Index<3> &indexCenter,
                                                                  std::vector<itk::Index<3>> &tumorSafetyMarginPixels,
                                                                  mitk::Image::Pointer &image,
                                                                  double &radius,
                                                                  mitk::Vector3D &ImageDimension,
                                                                  mitk::Vector3D &ImageSpacing);

  static double GetPercentageOfTumorvolumeInsideZone(itk::Index<3> &centerOfVolume,
                                                     mitk::Image::Pointer image,
                                                     mitk::Vector3D &imageSpacing,
                                                     mitk::Vector3D &imageDimension,
                                                     double radius);

  static double GetPercentageOfNonAblatedTumorvolumeInsideZone(itk::Index<3> &centerOfVolume,
                                                               mitk::Image::Pointer image,
                                                               mitk::Vector3D &imageSpacing,
                                                               mitk::Vector3D &imageDimension,
                                                               double radius);

  static int GetNumberOfAblatedPoints(itk::Index<3> &centerOfVolume,
                                      mitk::Image::Pointer image,
                                      mitk::Vector3D &imageSpacing,
                                      mitk::Vector3D &imageDimension,
                                      double radius);

  static std::vector<std::vector<itk::Index<3>>> FindAgglomerations(mitk::Image::Pointer image,
                                                                    mitk::Vector3D &imageSpacing,
                                                                    mitk::Vector3D &imageDimension);

  static std::vector<itk::Index<3>> FindFullAgglomeration(mitk::Image::Pointer image, itk::Index<3> &startingIndex);

  static bool CheckIfPixelIsElementOfAgglomeration(std::vector<itk::Index<3>> &agglomeration, itk::Index<3> index);

  static bool CheckIfPixelIsElementOfAgglomerationList(std::vector<std::vector<itk::Index<3>>> &agglomerationList,
                                                       itk::Index<3> pixel);

  static void SetSolutionValueStatistics(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans);

  static void SetMinMaxAblationZoneNumber(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans);

  static void SetMinMaxOverlapVolume(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans);

  static void SetMinMaxVolumeOutsideFactor(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans);

  // static void SetMinMaxAgglomerations(std::vector<mitk::AblationPlan::Pointer> AllFoundPlans);

  static int intRand(const int &min, const int &max);

private:
  AblationUtils();
  virtual ~AblationUtils();
};

#endif // ABLATIONUTILS_H
