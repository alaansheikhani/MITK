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
#include <mitkImage.h>

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

std::vector<std::vector<double>> RadiusModellingUtils::defineTimesAndPowers()
{
  std::vector<double> ablationTime;
  std::vector<double> ablationPower;
  std::vector<std::vector<double>> res;
  for (int i = 2; i < 37; i++)
  {
    ablationTime.push_back(60.0 + i * 15.0);
  }
  res.push_back(ablationTime);
  for (int i = 3; i < 11; i++)
  {
    ablationPower.push_back(30.0 + i * 5.0);
  }
  res.push_back(ablationPower);
  return res;
}

int RadiusModellingUtils::getMinLimitTime()
{
  std::vector<double> time = defineTimesAndPowers().at(0);
  double minTime = time.at(distance(time.begin(), min_element(time.begin(), time.end())));
  return minTime;
}

int RadiusModellingUtils::getMaxLimitTime()
{
  std::vector<double> time = defineTimesAndPowers().at(0);
  double maxTime = time.at(distance(time.begin(), max_element(time.begin(), time.end())));
  return maxTime;
}

int RadiusModellingUtils::getMinLimitPower()
{
  std::vector<double> power = defineTimesAndPowers().at(1);
  double minPower = power.at(distance(power.begin(), min_element(power.begin(), power.end())));
  return minPower;
}

int RadiusModellingUtils::getMaxLimitPower()
{
  std::vector<double> power = defineTimesAndPowers().at(1);
  double maxPower = power.at(distance(power.begin(), max_element(power.begin(), power.end())));
  return maxPower;
}

double RadiusModellingUtils::RDophi_Post(int time, int power)
{
  double R =
    4.312 + 0.007297 * time + 0.1 * power + 1.015 * pow(10, -5) * pow(time, 2) + 6.11 * pow(10, -19) * time * power;
  return R;
}

double RadiusModellingUtils::RDophi_Pre(int time, int power)
{
  double R = 3.339 - 0.004371 * time + 0.1864 * power + 1.961 * pow(10, -5) * pow(time, 2) + 0.0006476 * time * power -
             0.002331 * pow(power, 2) - 1.017 * pow(10, -7) * pow(time, 2) * power -
             4.389 * pow(10, -6) * time * pow(power, 2) + 2.218e-05 * pow(power, 3);
  return R;
}

double RadiusModellingUtils::DophiTimeSolved_Post(int power, double radius)
{
  // define the coefficients when solving the quadretic equation
  double a = 1.015 * pow(10, -5);
  double b = 0.007297 + 6.11 * pow(10, -19) * power;
  double c = 0.1 * power + 4.312 - radius;
  double time1, time2;
  // calculate the discriminant; the discriminant will always be positive (in range of the defined limits)
  double discriminant = b * b - 4 * a * c;
  // get the ablation time
  if (discriminant > 0)
  {
    time1 = (-b + sqrt(discriminant)) / (2 * a);
    time2 = (-b - sqrt(discriminant)) / (2 * a);
    // in the defined power and time limits of the model, time1 is always pos, whereas time2 is neg.
    if (time1 > 0 && time2 > 0)
    {
      if (time1 > time2)
      {
        return time1;
      }
      return time2;
    }
    else if (time1 > 0 || time2 > 0)
    {
      if (time1 > 0)
      {
        return time1;
      }
      return time2;
    }
    // MITK_INFO << "Ablation time cann't be calculated";
    return 0;
  }
  else if (discriminant == 0)
  {
    time1 = -b / (2 * a);
    return time1;
  }
  else
  {
    // MITK_INFO << "Ablation time cann't be calculated";
    return 0;
  }
}

double RadiusModellingUtils::DophiTimeSolved_Pre(int power, double radius)
{
  // define the coefficients when solving the quadretic equation
  double a = 1.961 * pow(10, -5) - 1.017 * pow(10, -7) * power;
  double b = 0.0006476 * power - 0.004371;
  double c = 3.339 + 0.1864 * power - 0.002331 * pow(power, 2) - 2.218 * pow(10, -5) * pow(power, 3) - radius;
  double time1, time2;
  // calculate the discriminant; the discriminant will always be positive (in range of the defined limits)
  double discriminant = b * b - 4 * a * c;
  // get the ablation time
  if (discriminant > 0)
  {
    time1 = (-b + sqrt(discriminant)) / (2 * a);
    time2 = (-b - sqrt(discriminant)) / (2 * a);
    // in the defined power and time limits of the model, time1 is always pos, whereas time2 is neg.
    if (time1 > 0 && time2 > 0)
    {
      if (time1 > time2)
      {
        return time1;
      }
      return time2;
    }
    else if (time1 > 0 || time2 > 0)
    {
      if (time1 > 0)
      {
        return time1;
      }
      return time2;
    }
    // MITK_INFO << "Ablation time cann't be calculated";
    return 0;
  }
  else if (discriminant == 0)
  {
    time1 = -b / (2 * a);
    return time1;
  }
  else
  {
    // MITK_INFO << "Ablation time cann't be calculated";
    return 0;
  }
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
  double minS = shrinkage(getMinLimitTime(), getMinLimitPower());
  return minS;
}

double RadiusModellingUtils::getMaxShrinkage()
{
  double maxS = shrinkage(getMaxLimitTime(), getMaxLimitPower());
  return maxS;
}

double RadiusModellingUtils::getMinRadiusDophi()
{
  double minR = RDophi_Post(getMinLimitTime(), getMinLimitPower());
  return minR;
}

double RadiusModellingUtils::getMaxRadiusDophi()
{
  double maxR = RDophi_Post(getMaxLimitTime(), getMaxLimitPower());
  return maxR;
}

double RadiusModellingUtils::getMinRadiusEmprint()
{
  double minR = REmprint(getMinLimitTime(), getMinLimitPower());
  return minR;
}

double RadiusModellingUtils::getMaxRadiusEmprint()
{
  double maxR = REmprint(getMaxLimitTime(), getMaxLimitPower());
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

std::vector<double> RadiusModellingUtils::calculateTimeAndPowreOfPostRadiusDophi(double radius)
{
  std::vector<std::vector<double>> availableSettings = defineTimesAndPowers();
  double ablationPower = availableSettings.at(1).at(3);
  double ablationTime = DophiTimeSolved_Post(ablationPower, radius);
  std::vector<double> ablationSettings;
  if (ablationTime >= getMinLimitTime() && ablationTime <= getMaxLimitTime())
  {
    ablationSettings.push_back(ablationTime);
    ablationSettings.push_back(ablationPower);
  }
  else if (ablationTime < getMinLimitTime())
  {
    for (int i = 1; i < size(availableSettings.at(1)) - 3; i++)
    {
      ablationPower = availableSettings.at(1).at(6) - 5.0 * i;
      ablationTime = DophiTimeSolved_Post(ablationPower, radius);
      if (ablationTime >= getMinLimitTime() && ablationTime <= getMaxLimitTime())
      {
        ablationSettings.push_back(ablationTime);
        ablationSettings.push_back(ablationPower);
      }
    }
  }
  else
  {
    for (int i = 1; i < size(availableSettings.at(1)) - 5; i++)
    {
      ablationPower = availableSettings.at(1).at(6) + 5.0 * i;
      ablationTime = DophiTimeSolved_Post(ablationPower, radius);
      if (ablationTime >= getMinLimitTime() && ablationTime <= getMaxLimitTime())
      {
        ablationSettings.push_back(ablationTime);
        ablationSettings.push_back(ablationPower);
      }
    }
  }
  return ablationSettings;
}

std::vector<double> RadiusModellingUtils::calculateTimeAndPowreOfPreRadiusDophi(double radius)
{
  std::vector<std::vector<double>> availableSettings = defineTimesAndPowers();
  double ablationPower = availableSettings.at(1).at(3);
  double ablationTime = DophiTimeSolved_Pre(ablationPower, radius);
  std::vector<double> ablationSettings;
  if (ablationTime >= getMinLimitTime() && ablationTime <= getMaxLimitTime())
  {
    ablationSettings.push_back(ablationTime);
    ablationSettings.push_back(ablationPower);
  }
  else if (ablationTime < getMinLimitTime())
  {
    for (int i = 1; i < size(availableSettings.at(1)) - 3; i++)
    {
      ablationPower = availableSettings.at(1).at(6) - 5.0 * i;
      ablationTime = DophiTimeSolved_Pre(ablationPower, radius);
      if (ablationTime >= getMinLimitTime() && ablationTime <= getMaxLimitTime())
      {
        ablationSettings.push_back(ablationTime);
        ablationSettings.push_back(ablationPower);
      }
    }
  }
  else
  {
    for (int i = 1; i < size(availableSettings.at(1)) - 5; i++)
    {
      ablationPower = availableSettings.at(1).at(6) + 5.0 * i;
      ablationTime = DophiTimeSolved_Pre(ablationPower, radius);
      if (ablationTime >= getMinLimitTime() && ablationTime <= getMaxLimitTime())
      {
        ablationSettings.push_back(ablationTime);
        ablationSettings.push_back(ablationPower);
      }
    }
  }
  return ablationSettings;
}

double RadiusModellingUtils::calculateShrinkageOfPostRadiusDophi(double radius)
{
  std::vector<double> ablationSettings = calculateTimeAndPowreOfPostRadiusDophi(radius);
  double shrinkageValue = shrinkage(ablationSettings.at(0), ablationSettings.at(1));
  return shrinkageValue;
}

double RadiusModellingUtils::calculateShrinkageOfPreRadiusDophi(double radius)
{
  std::vector<double> ablationSettings = calculateTimeAndPowreOfPreRadiusDophi(radius);
  double shrinkageValue = shrinkage(ablationSettings.at(0), ablationSettings.at(1));
  return shrinkageValue;
}
