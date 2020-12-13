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

#ifndef ABLATIONPLANOPTIMIZER_H
#define ABLATIONPLANOPTIMIZER_H
#include "mitkAblationPlan.h"
#include "mitkImagePixelReadAccessor.h"
#include "mitkImagePixelWriteAccessor.h"

/**
  \brief TODO

  \ingroup ${plugin_target}_internal
*/
namespace mitk
{
  class AblationPlanOptimizer
  {
  public:
    static void RemoveNotNeededVolumes(mitk::AblationPlan::Pointer plan);
    //static void ShrinkAndMoveZones(mitk::AblationPlan::Pointer plan);
    static void Optimize(mitk::AblationPlan::Pointer plan, std::vector<mitk::AblationZone> &tempAblationZones);
    static void MoveZonesTowardsCenter(mitk::AblationPlan::Pointer plan,
                                       std::vector<mitk::AblationZone> &tempAblationZones);
    static void CheckIfDistanceBeetweenZonesTooSmall(mitk::AblationPlan::Pointer plan,
                                                     std::vector<mitk::AblationZone> &tempAblationZones);
    //static bool CheckNewCenterForNonAblatedVolume(mitk::AblationPlan::Pointer plan,
    //                                              int zoneNumber,
    //                                              itk::Index<3> newCenter);
    //static void SetNewAblationVolume(mitk::AblationPlan::Pointer plan, int zoneNumber, itk::Index<3> newCenter);

  private:
    AblationPlanOptimizer();
    virtual ~AblationPlanOptimizer();
  };
} // namespace mitk

#endif // ABLATIONPLANOPTIMIZER_H
