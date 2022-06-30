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

#include "mitkAblationPlanningLogging.h"
#include <mitkSceneIO.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <QString>

mitk::AblationPlanningLogging::AblationPlanningLogging(){
}

mitk::AblationPlanningLogging::~AblationPlanningLogging(){

}

void mitk::AblationPlanningLogging::WriteHeader(){
  std::ofstream file;
  file.open(m_FileName.c_str(), std::ios_base::app); // append to file
  file
    << "Case Name; Tumor Volume [ml]; Safety Margin [mm]; Tumor + Safetymargin Volume [ml]; "          // 1,2,3,4
    << "Param: Minimal Radius [mm]; Param: Desired Radius [mm]; Param: Maximal Radius [mm]; "          // 5,6,7
    << "Param: Shrinking [percentage]; Param: Tolerance Non-Ablated Volume [percentage]; Iterations; " // 8,9,10
    << "Number of Zones; Factor Non-Ablated Volume; Factor Overlapping Zones; Factor Ablated Volume Outside Of Tumor + "
       "Safetymargin Volume; " // 11,12,13,14
    << "Total Ablation Volume; Solution Value from Metric; Radii of all Zones; Radii of all zones after shrinking;" // 15,16,17,18
    << "Distances from tumor COG to zones centers [mm]; Tumor COG in world coordinates; Position of zone centers "
       "relative to COG\n"; // 19,20,21
  file.close();
}

void mitk::AblationPlanningLogging::WriteDataSet(
  mitk::AblationPlan::Pointer plan,
  mitk::DataNode::Pointer tumorNode,
  mitk::AblationPlanningLogging::AblationPlanningParameterSet parameterSet,
  std::string name)
{
  std::ofstream file;
  file.open(m_FileName.c_str(), std::ios_base::app);                                // append to file
  file << name << ";"                                                               // 1
       << plan->GetStatistics().tumorVolume << ";"                                  // 2
       << parameterSet.safetyMargin << ";"                                          // 3
       << plan->GetStatistics().tumorAndSafetyMarginVolume << ";"                   // 4
       << parameterSet.minRadius << ";"                                             // 5
       << parameterSet.desiredRadius << ";"                                         // 6
       << parameterSet.maxRadius << ";"                                             // 7
       << parameterSet.tissueShrinking << ";"                                       // 8
       << parameterSet.toleranceNonAblatedVolume << ";"                             // 9
       << parameterSet.iterations << ";"                                            // 10
       << plan->GetNumberOfZones() << ";"                                           // 11
       << plan->GetStatistics().factorNonAblatedVolume / 100 << ";"                 // 12
       << plan->GetStatistics().factorOverlappingAblationZones / 100 << ";"         // 13
       << plan->GetStatistics().factorAblatedVolumeOutsideSafetyMargin / 100 << ";" // 14
       << plan->GetStatistics().totalAblationVolume << ";"                          // 15
       << plan->GetSolutionValue() << ";";                                          // 16
  for (int j = 0; j < plan->GetNumberOfZones(); j++)
  {
    file << plan->GetAblationZone(j)->radius << "-"; // 17
  }
  file << ";";
  for (int j = 0; j < plan->GetNumberOfZones(); j++)
  {
    file << plan->GetAblationZone(j)->radius / (1 + ((double)parameterSet.tissueShrinking / 100.0)) << "-"; // 18
  }
  file << ";";
  distances = AblationUtils::CalculateDistancesOfCOGToCenters(parameterSet.COGravity, plan);
  for (int j = 0; j < distances.size(); j++)
  {
    file << distances[j] << "-"; // 19
  }
  file << ";"
       << "[" << parameterSet.COGravity->GetPoint(0)[0] << "," << parameterSet.COGravity->GetPoint(0)[1] << "," // 20
       << parameterSet.COGravity->GetPoint(0)[2] << "]"
       << ";";
  for (int index = 1; index < parameterSet.relativeCoordinates->GetSize(); ++index)
  {
    file << "[" << parameterSet.relativeCoordinates->GetPoint(index)[0] << ","
         << parameterSet.relativeCoordinates->GetPoint(index)[1] << "," // 21
         << parameterSet.relativeCoordinates->GetPoint(index)[2] << "] | ";
  }
  file << "\n";
  file.close();
}

void mitk::AblationPlanningLogging::WriteScene(mitk::DataStorage::Pointer storage, std::string name)
{
  mitk::SceneIO::Pointer mySceneIO = mitk::SceneIO::New();
  QString filenameScene = QString(m_FileName.c_str()) + QString(name.c_str()) + "_mitkScene.mitk";
  mitk::NodePredicateNot::Pointer isNotHelperObject =
    mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true)));
  mitk::DataStorage::SetOfObjects::ConstPointer nodesToBeSaved = storage->GetSubset(isNotHelperObject);
  mySceneIO->SaveScene(nodesToBeSaved, storage, filenameScene.toStdString().c_str());
}
