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

#ifndef MITKABLATIONPLANNINGALGORITHM_H_HEADER_INCLUDED_
#define MITKABLATIONPLANNINGALGORITHM_H_HEADER_INCLUDED_
#include <itkObject.h>
#include <mitkCommon.h>
#include <mitkImage.h>
#include "mitkAblationUtils.h"

namespace mitk
{
  /**Documentation
   * \brief AblationPlanningAlgorithm
   */
  class AblationPlanningAlgorithm : public itk::Object
  {
  public:
    mitkClassMacroItkParent(AblationPlanningAlgorithm, itk::Object);
    itkFactorylessNewMacro(Self);

  protected:
    AblationPlanningAlgorithm();

    ~AblationPlanningAlgorithm() override;

  private:

  };
} // namespace mitk
#endif /* MITKABLATIONPLANNINGALGORITHM_H_HEADER_INCLUDED_ */
