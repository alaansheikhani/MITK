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
#include "mitkAblationPlan.h"
#include <mitkDataNode.h>
#include <mitkDataStorage.h>
#include <vector>

namespace mitk
{
  /**Documentation
   * \brief AblationPlanLogging
   */
  class AblationPlanningLogging : public itk::Object
  {
  public:
    struct AblationPlanningParameterSet{
      double minRadius;
      double desiredRadius;
      double maxRadius;
      double tissueShrinking;
      double toleranceNonAblatedVolume;
      double safetyMargin;
      int iterations;
    };

    mitkClassMacroItkParent(AblationPlanningLogging, itk::DataObject);
    itkFactorylessNewMacro(Self);
    itkSetMacro(FileName,std::string);
    void WriteHeader();
    /** Writes one line with information about that dataset to the logging file. Don't forget
        to write an heade before if the file is empty before.*/
    void WriteDataSet(mitk::AblationPlan::Pointer plan, mitk::DataNode::Pointer tumorNode, mitk::AblationPlanningLogging::AblationPlanningParameterSet parameterSet, std::string name);
    /** Writs a whole planning scene to a subfolder with the name of this
        data set. */
    void WriteScene(mitk::DataStorage::Pointer storage, std::string name);


  protected:
    AblationPlanningLogging();

    ~AblationPlanningLogging() override;

  private:
    std::string m_FileName;

  };
} // namespace mitk
#endif /* MITKABLATIONPLANLOG_H_HEADER_INCLUDED_ */
