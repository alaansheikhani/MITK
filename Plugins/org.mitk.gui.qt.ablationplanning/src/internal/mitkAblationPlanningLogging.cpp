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

mitk::AblationPlanningLogging::AblationPlanningLogging(){
}

mitk::AblationPlanningLogging::~AblationPlanningLogging(){

}

void mitk::AblationPlanningLogging::WriteHeader(){
  std::ofstream file;
  file.open(m_FileName.c_str(), std::ios_base::app); //append to file
  file << "Tumor Name,Tumor Volume,Tumor + Safetymargin Volume,Plan Nr,NumberOfZones,Factor Non-Ablated Volume,Factor "
            "Overlapping Zones,Factor Ablated Volume Outside Of Tumor + Safetymargin Volume,Total Ablation "
            "Volume,SolutionValue,Radius "
      "Of Zones\n";
  file.close();
}

void mitk::AblationPlanningLogging::WriteDataSet(mitk::AblationPlan::Pointer plan, mitk::DataNode::Pointer tumorNode, std::string name){
  std::ofstream file;
  file.open(m_FileName.c_str(), std::ios_base::app); //append to file
  file   << "Name" << ","
         << "Volume" << ","
         << "VolumeWSM" << ","
         << "Final" << ","
         << plan->GetNumberOfZones() << ","
         << plan->GetStatistics().factorNonAblatedVolume / 100 << ","
         << plan->GetStatistics().factorOverlappingAblationZones / 100 << ","
         << plan->GetStatistics().factorAblatedVolumeOutsideSafetyMargin / 100 << ","
         << plan->GetStatistics().totalAblationVolume << "," << plan->GetSolutionValue()
         << ",";
  for (int j = 0; j < plan->GetNumberOfZones(); j++){
        file << plan->GetAblationZone(j)->radius << " ";
  }
  file << "\n";
  file.close();
}

void mitk::AblationPlanningLogging::WriteScene(mitk::DataStorage::Pointer storage, std::string name){

}
