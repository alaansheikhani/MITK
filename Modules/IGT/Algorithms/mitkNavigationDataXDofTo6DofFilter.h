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

#ifndef NAVIGATIONDATAXDOFTO6DOFFILTER_H
#define NAVIGATIONDATAXDOFTO6DOFFILTER_H

#include "MitkIGTExports.h"
#include "mitkNavigationDataToNavigationDataFilter.h"

namespace mitk
{
 /**Documentation
  * \brief NavigationDataXDofTo6DofFilter merges xDoF (e.g. 5DoF) sensors to a 6DoF sensor, while emits multiple mitk::Message messages when the input NavigationData values
  * change
  *
  * This filter can have multiple inputs. It emits
  * the following Messages if an input navigation data values changed since the last Update()
  * - PositionChangedMessage
  * - OrientationChangedMessage
  * - ErrorChangedMessage
  * - TimeStampChangedMessage
  * - DataValidChangedMessage
  *
  * The first parameter of these messages is the new value, the second is the index of the input that has changed
  * The filter has as many outputs as it has inputs. It copies the inputs to the outputs after sending the messages.
  * \ingroup IGT
  */
  class MITKIGT_EXPORT NavigationDataXDofTo6DofFilter : public NavigationDataToNavigationDataFilter
  {
  public:
    struct Landmark
    {
      mitk::Point3D source_pt;
      mitk::Point3D target_pt_sensor_coordinates;
      unsigned int inputID;
      unsigned int outputID;
    };

    mitkClassMacro(NavigationDataXDofTo6DofFilter, NavigationDataToNavigationDataFilter);
    itkNewMacro(Self);

    // void InitializeTransform(int inputID, int outputID, mitk::Point3D originPoint);
    void AddLandmarkFor6DoF(mitk::Point3D source_pt,
                    mitk::Point3D target_pt,
                    unsigned int inputID,
                    unsigned int outputID); // just one sourcePoint for one Input (origin -> (0,0,0))

  protected:
    NavigationDataXDofTo6DofFilter();
    ~NavigationDataXDofTo6DofFilter() override;
    void GenerateData() override;
    int GetNumberOfSourcePointsForOutput(unsigned int outputID);
    void PassThrough(unsigned int inputID, unsigned int outputID);
    void SetLandmarksForOutput(unsigned int outputID);

    std::vector<Landmark> m_Landmarks;
    vtkSmartPointer<vtkPoints> m_SourcePoints;
    vtkSmartPointer<vtkPoints> m_TargetPoints;

  private:
    std::vector<int> m_OutputList; ///> Holds all output ids from landmarks once
    void UpdateListOfAllOutputs(); ///> Updates the output id list from the list of landmarks
    bool IsOutputInList(int outputID); ///> Checks if an output id is already in the list

  };

} // namespace mitk




#endif // NAVIGATIONDATAXDOFTO6DOFFILTER_H
