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
  static double RDophi(int time, int power);
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

private:
  RadiusModellingUtils();
  virtual ~RadiusModellingUtils();
};

#endif
