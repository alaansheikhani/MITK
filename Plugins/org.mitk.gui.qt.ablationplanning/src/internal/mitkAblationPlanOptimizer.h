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
      static void Optimize(mitk::AblationPlan::Pointer plan, std::vector<mitk::AblationZone> &tempAblationZones);

    private:
      AblationPlanOptimizer();
      virtual ~AblationPlanOptimizer();
};
}

#endif // ABLATIONPLANOPTIMIZER_H
