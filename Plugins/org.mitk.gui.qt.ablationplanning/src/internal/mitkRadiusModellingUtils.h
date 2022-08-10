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

#ifndef RADIUSMODELLINGUTILS_H
#define RADIUSMODELLINGUTILS_H

/**
  \brief RadiusModellingUtils static class for doing radius calculation of the min and max radius of selected MWA System

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \ingroup ${plugin_target}_internal
*/

class RadiusModellingUtils
{
public:
  static std::vector<std::vector<double>> defineTimesAndPowers();
  static int getMinLimitTime();
  static int getMaxLimitTime();
  static int getMinLimitPower();
  static int getMaxLimitPower();
  static double RDophi_Post(int time, int power);
  static double RDophi_Pre(int time, int power);
  static double DophiTimeSolved_Post(int power, double radius);
  static double DophiTimeSolved_Pre(int power, double radius);
  static double REmprint(int time, int power);
  static double shrinkage(int time, int power);
  static double getMinShrinkage();
  static double getMaxShrinkage();
  static double getMinRadiusDophi();
  static double getMaxRadiusDophi();
  static double getMinRadiusEmprint();
  static double getMaxRadiusEmprint();
  static double getPreAblationMinRadiusDophi();
  static double getPreAblationMaxRadiusDophi();
  static double getPreAblationMinRadiusEmprint();
  static double getPreAblationMaxRadiusEmprint();
  static std::vector<double> calculateTimeAndPowreOfPostRadiusDophi(double radius);
  static std::vector<double> calculateTimeAndPowreOfPreRadiusDophi(double radius);
  static double calculateShrinkageOfPostRadiusDophi(double radius);
  static double calculateShrinkageOfPreRadiusDophi(double radius);

private:
  RadiusModellingUtils();
  virtual ~RadiusModellingUtils();
};

#endif
