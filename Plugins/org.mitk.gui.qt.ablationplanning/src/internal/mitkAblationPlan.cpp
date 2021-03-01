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

#include "mitkAblationPlan.h"

mitk::AblationPlan::AblationPlan()
  : m_SegmentationImage(NULL),
    m_AblationZones(std::vector<mitk::AblationZone>()),
    m_StatsSet(false),
    m_Stats(mitk::AblationPlan::AblationPlanStatistics())
{
}

mitk::AblationPlan::~AblationPlan() {}

mitk::AblationPlan::AblationPlanStatistics mitk::AblationPlan::GetStatistics()
{
  return m_Stats;
}

void mitk::AblationPlan::SetStatistics(mitk::AblationPlan::AblationPlanStatistics s)
{
  m_StatsSet = true;
  m_Stats = s;
}

bool mitk::AblationPlan::AddAblationZone(mitk::AblationZone newZone)
{
  m_AblationZones.push_back(newZone);
  return true;
}

int mitk::AblationPlan::GetNumberOfZones()
{
  return m_AblationZones.size();
}

mitk::AblationZone *mitk::AblationPlan::GetAblationZone(int id)
{
  return &m_AblationZones.at(id);
}

bool mitk::AblationPlan::RemoveAblationZone(mitk::AblationZone &zone)
{
  for (int i = 0; i < GetNumberOfZones(); i++)
  {
    if (m_AblationZones.at(i).IsEqual(zone))
    {
      m_AblationZones.erase(m_AblationZones.begin() + i);
      return true;
    }
  }
  return false;
}

bool mitk::AblationPlan::RemoveAblationZone(int id)
{
  m_AblationZones.erase(m_AblationZones.begin() + id);
  return true;
}

void mitk::AblationPlan::SetSolutionValue(double newValue)
{
  this->m_SolutionValue = newValue;
}

double mitk::AblationPlan::GetSolutionValue()
{
  return this->m_SolutionValue;
}

//----------- comparison: compare m_solutionValue
int mitk::AblationPlan::CompareTo(mitk::AblationPlan::Pointer b)
{
  if (this->GetSolutionValue() > b->GetSolutionValue())
  {
    MITK_INFO << "Improved from " << b->GetNumberOfZones() << " zones to " << this->GetNumberOfZones() << " zones.";
    return -1;
  }
  else if (this->GetSolutionValue() < b->GetSolutionValue())
  {
    return 1;
  }
  else if (this->GetSolutionValue() == b->GetSolutionValue())
  {
    MITK_INFO << "Solutions are Equal";
    return 0;
  }
  return 0;
}

double mitk::AblationPlan::CalculateSolutionValueOfZoneNumbers() // 0 worst 1 best
{
  return (double(this->GetMaxAblationZoneNumber()) - double(this->GetNumberOfZones())) /
         (double(this->GetMaxAblationZoneNumber()) - double(this->GetMinAblationZoneNumber()));
}

double mitk::AblationPlan::CalculateSolutionValueOfZoneDifference()
{
  std::vector<int> zoneRadius{};
  int counter{0};
  for (size_t i = 0; i < this->GetNumberOfZones(); i++)
  {
    zoneRadius.push_back(this->GetAblationZone(i)->radius);
  }
  std::sort(zoneRadius.begin(), zoneRadius.end(), std::greater<>());
  for (int i = 1; i < zoneRadius.size(); i++)
  {
    if (zoneRadius.at(i - 1) != zoneRadius.at(i))
    {
      counter++;
    }
  }
  return double(counter) / double(zoneRadius.size() - 1);
}

double mitk::AblationPlan::CalculateSolutionValueOfOverlappingZones()
{
  return (this->GetMaxOverlappingZonesFactor() - this->GetStatistics().factorOverlappingAblationZones) /
         (this->GetMaxOverlappingZonesFactor() - this->GetMinOverlappingZonesFactor());
}

double mitk::AblationPlan::CalculateSolutionValueOfVolumeOutsideFactor()
{
  return (this->GetMaxVolumeOutsideFactor() - this->GetStatistics().factorAblatedVolumeOutsideSafetyMargin) /
         (this->GetMaxOverlappingZonesFactor() - this->GetMinOverlappingZonesFactor());
}

double mitk::AblationPlan::CalculateSolutionValueOfNonAblated()
{
  return ((m_Stats.factorMaxNonAblatedVolume * 100) - this->GetStatistics().factorNonAblatedVolume) /
         (m_Stats.factorMaxNonAblatedVolume * 100);
}

void mitk::AblationPlan::CalculcateSolutionValue()
{
  double valueNumberOfZones = this->CalculateSolutionValueOfZoneNumbers();
  double valueNonAblated = this->CalculateSolutionValueOfNonAblated();
  double valueRadius = this->CalculateSolutionValueOfZoneDifference();
  double valueOverlap = this->CalculateSolutionValueOfOverlappingZones();
  double valueHealthy = this->CalculateSolutionValueOfVolumeOutsideFactor();

  double solutionValue = (valueNumberOfZones * 50 + valueHealthy * 50) / 100;
  if (this->GetStatistics().factorNonAblatedVolume > this->m_Stats.factorMaxNonAblatedVolume * 100)
  {
    solutionValue -= 1000.0;
  }

  SetSolutionValue(solutionValue);
}

void mitk::AblationPlan::SetMaxAblationZoneNumber(int n)
{
  this->m_MaxAblationZoneNumber = n;
}

void mitk::AblationPlan::SetMinAblationZoneNumber(int n)
{
  this->m_MinAblationZoneNumber = n;
}

int mitk::AblationPlan::GetMaxAblationZoneNumber()
{
  return this->m_MaxAblationZoneNumber;
}

int mitk::AblationPlan::GetMinAblationZoneNumber()
{
  return this->m_MinAblationZoneNumber;
}

void mitk::AblationPlan::SetMaxOverlappingZonesFactor(double overlapFactor)
{
  this->m_MaxOverlappingZonesFactor = overlapFactor;
}
void mitk::AblationPlan::SetMinOverlappingZonesFactor(double overlapFactor)
{
  this->m_MinOverlappingZonesFactor = overlapFactor;
}

double mitk::AblationPlan::GetMaxOverlappingZonesFactor()
{
  return this->m_MaxOverlappingZonesFactor;
}

double mitk::AblationPlan::GetMinOverlappingZonesFactor()
{
  return this->m_MinOverlappingZonesFactor;
}

void mitk::AblationPlan::SetSegmentationImage(mitk::Image::Pointer s)
{
  this->m_SegmentationImage = s;
}

mitk::Image::Pointer mitk::AblationPlan::GetSegmentationImage()
{
  return m_SegmentationImage;
}

void mitk::AblationPlan::SetMaxVolumeOutsideFactor(double factorOutside)
{
  this->m_MaxVolumeOutsideFactor = factorOutside;
}
void mitk::AblationPlan::SetMinVolumeOutsideFactor(double factorOutside)
{
  this->m_MinVolumeOutsideFactor = factorOutside;
}
double mitk::AblationPlan::GetMaxVolumeOutsideFactor()
{
  return this->m_MaxVolumeOutsideFactor;
}
double mitk::AblationPlan::GetMinVolumeOutsideFactor()
{
  return this->m_MinVolumeOutsideFactor;
}