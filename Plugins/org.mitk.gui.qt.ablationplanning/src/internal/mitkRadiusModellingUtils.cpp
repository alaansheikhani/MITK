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
#include "mitkRadiusModellingUtils.h"

// Qt
#include <QMessageBox>

// rng

#if defined(_MSC_VER) // Visual studio
#define thread_local __declspec(thread)
#elif defined(__GCC__) // GCC
#define thread_local __thread
#endif

//=====================Konstruktor/Destruktor===================================
RadiusModellingUtils::RadiusModellingUtils() {}

RadiusModellingUtils::~RadiusModellingUtils() {}

double RadiusModellingUtils::RDophi(int time, int power)
{
  double R =
    4.312 + 0.007297 * time + 0.1 * power + 1.015 * pow(10, -5) * pow(time, 2) + 6.11 * pow(10, -19) * time * power;
  return R;
}

double RadiusModellingUtils::REmprint(int time, int power)
{
  double R = -0.006234 + 0.05341 * time + 0.05382 * power - 6.961 * pow(10, -5) * pow(time, 2) -
             1.718 * pow(10, -5) * time * power + 3.961 * pow(10, -8) * pow(time, 3) + 2.718 * pow(time, 2) * power;
  return R;
}

double RadiusModellingUtils::shrinkage(int time, int power)
{
  double S = ((0.1536 * power + 56.02) /
              (1 + pow((time / (611.7 + pow(150.7, 2) * exp(-0.1491 * power))), (0.01165 * power - 1.374)))) -
             0.476;
  return S / 100;
}

double RadiusModellingUtils::getMinShrinkage()
{
  double minS = shrinkage(60, 20);
  return minS;
}

double RadiusModellingUtils::getMaxShrinkage()
{
  double maxS = shrinkage(600, 80);
  return maxS;
}

double RadiusModellingUtils::getMinRadiusDophi()
{
  double minR = RDophi(60, 20);
  return minR;
}

double RadiusModellingUtils::getMaxRadiusDophi()
{
  double maxR = RDophi(600, 80);
  return maxR;
}

double RadiusModellingUtils::getMinRadiusEmprint()
{
  double minR = REmprint(60, 20);
  return minR;
}

double RadiusModellingUtils::getMaxRadiusEmprint()
{
  double maxR = REmprint(600, 80);
  return maxR;
}

double RadiusModellingUtils::getPreAblationMinRadiusDophi()
{
  double RPre = getMinRadiusDophi() / (1 - getMinShrinkage());
  return RPre;
}

double RadiusModellingUtils::getPreAblationMaxRadiusDophi()
{
  double RPre = getMaxRadiusDophi() / (1 - getMaxShrinkage());
  return RPre;
}

double RadiusModellingUtils::getPreAblationMinRadiusEmprint()
{
  double RPre = getMinRadiusEmprint() / (1 - getMinShrinkage());
  return RPre;
}

double RadiusModellingUtils::getPreAblationMaxRadiusEmprint()
{
  double RPre = getMaxRadiusEmprint() / (1 - getMaxShrinkage());
  return RPre;
}