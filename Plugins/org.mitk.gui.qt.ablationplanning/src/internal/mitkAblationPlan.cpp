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
  : m_SegmentationImage(NULL), m_AblationZones(std::vector<mitk::AblationZone>())
{
}

mitk::AblationPlan::~AblationPlan() {}

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

bool mitk::AblationPlan::RemoveAblationZone(int id){
  m_AblationZones.erase(m_AblationZones.begin() + id);
  return true;
}

int mitk::AblationPlan::CompareTo(mitk::AblationPlan::Pointer b)
{
  if (this->GetNumberOfZones() > b->GetNumberOfZones())
  {
    return 1;
  }
  else if (this->GetNumberOfZones() < b->GetNumberOfZones())
  {
    return -1;
  }
  else
  {
    //todo: other criteria!
    return 0;
  }
}
