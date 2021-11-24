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

#ifndef MITKABLATIONPLANLOG_H_HEADER_INCLUDED_
#define MITKABLATIONPLANLOG_H_HEADER_INCLUDED_
#include <mitkCommon.h>
#include <mitkImage.h>
#include "mitkAblationZone.h"
#include <vector>

namespace mitk
{
  /**Documentation
   * \brief AblationPlanLogging
   */
  class AblationPlanningLogging : public itk::Object
  {
  public:
    mitkClassMacroItkParent(AblationPlanningLogging, itk::DataObject);
    itkFactorylessNewMacro(Self);

  protected:
    AblationPlanningLogging();

    ~AblationPlanningLogging() override;

  private:

  };
} // namespace mitk
#endif /* MITKABLATIONPLANLOG_H_HEADER_INCLUDED_ */
