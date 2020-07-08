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

#include "mitkNavigationDataXDofTo6DofFilter.h"
#include <vtkLandmarkTransform.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkPoints.h>

mitk::NavigationDataXDofTo6DofFilter::NavigationDataXDofTo6DofFilter() : m_OutputList(std::vector<int>()), landmarks(std::vector<Landmark>()){
}

mitk::NavigationDataXDofTo6DofFilter::~NavigationDataXDofTo6DofFilter() {
}


void mitk::NavigationDataXDofTo6DofFilter::GenerateData()
{
  MITK_INFO << "Generate Data";
  //update list of all 6 DoF outputs
  UpdateListOfAllOutputs();

  //now compute all 6 DoF outputs
  for (int output_id : m_OutputList)
  {
    if (this->GetNumberOfSourcePointsForOutput(output_id) < 3)
    {
      MITK_WARN << "Cannot compute 6 DoF for output id " << output_id << ": SourcePointSet must contain at least 3 points";
      break;
    }
    else
    {
      this->SetLandmarksForOutput(output_id);

      // Compute transform
      vtkSmartPointer<vtkLandmarkTransform> transform = vtkSmartPointer<vtkLandmarkTransform>::New();
      transform->SetSourceLandmarks(m_SourcePoints);
      transform->SetTargetLandmarks(m_TargetPoints);
      transform->SetModeToRigidBody();
      transform->Modified();
      transform->Update();

      // Convert from vtk to itk data types
      itk::Matrix<float, 3, 3> rotationFloat = itk::Matrix<float, 3, 3>();
      itk::Vector<float, 3> translationFloat = itk::Vector<float, 3>();
      itk::Matrix<double, 3, 3> rotationDouble = itk::Matrix<double, 3, 3>();
      itk::Vector<double, 3> translationDouble = itk::Vector<double, 3>();
      vtkSmartPointer<vtkMatrix4x4> m = transform->GetMatrix();
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
        {
          rotationFloat[k][l] = m->GetElement(k, l);
          rotationDouble[k][l] = m->GetElement(k, l);
        }
      for (int k = 0; k < 3; k++)
      {
        translationFloat[k] = m->GetElement(k, 3);
        translationDouble[k] = m->GetElement(k, 3);
      }
      // Create affine transform 3D surface
      mitk::AffineTransform3D::Pointer mitkTransform = mitk::AffineTransform3D::New();
      mitkTransform->SetMatrix(rotationDouble);
      mitkTransform->SetOffset(translationDouble);

      //############### object is transformed ##########################
      // Save transform result
      mitk::NavigationData::Pointer transformedND = mitk::NavigationData::New(mitkTransform);
      transformedND->SetName("Computed 6DoF Output");
      mitk::NavigationData *output = this->GetOutput(output_id);
      output->Graft(transformedND); // copy all information from result to output
    }
  }

  //finally: pass through all other inputs
  for (unsigned int id=0; id<this->GetNumberOfInputs(); id++){
    if (!IsOutputInList(id)){
      PassThrough(id, id);
    }
  }
}

int mitk::NavigationDataXDofTo6DofFilter::GetNumberOfSourcePointsForOutput(unsigned int outputID)
{
  int count = 0;
  for (int i = 0; i < landmarks.size(); i++)
  {

    Landmark landmark = landmarks.at(i);
    if (landmark.outputID == outputID)
    {
      count++;
    }

  }
  return count;
}

void mitk::NavigationDataXDofTo6DofFilter::SetLandmarksForOutput(unsigned int outputID)
{
  for (int i = 0; i < landmarks.size(); i++)
  {
    vtkSmartPointer<vtkPoints> sourcePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> targetPoints = vtkSmartPointer<vtkPoints>::New();

    Landmark landmark = landmarks.at(i);
    if (landmark.outputID == outputID)
    {
      mitk::Point3D srcPt = landmark.source_pt;
      sourcePoints->InsertNextPoint(srcPt[0], srcPt[1], srcPt[2]);

      mitk::NavigationData::Pointer nd = mitk::NavigationData::New();
      nd->Graft(this->GetInput(landmark.inputID));


      mitk::NavigationData::Pointer targetPoint = mitk::NavigationData::New();
      targetPoint->SetPosition(landmark.target_pt_sensor_coordinates);
      targetPoint->Compose(nd);

      mitk::Point3D trgPt = targetPoint->GetPosition();
      targetPoints->InsertNextPoint(trgPt[0], trgPt[1], trgPt[2]);

    }
    m_SourcePoints = sourcePoints;
    m_TargetPoints = targetPoints;

  }
}

void mitk::NavigationDataXDofTo6DofFilter::AddLandmarkFor6DoF(mitk::Point3D source_pt,
                                                      mitk::Point3D target_pt_sensor_coordinates,
                                                      unsigned int inputID,
                                                      unsigned int outputID)
{
  vtkSmartPointer<vtkPoints> sourcePoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> targetPoints = vtkSmartPointer<vtkPoints>::New();

  Landmark landmark;
  landmark.inputID = inputID;
  landmark.outputID = outputID;
  landmark.source_pt = source_pt;
  landmark.target_pt_sensor_coordinates = target_pt_sensor_coordinates;

  // inserts the new landmark instance at front
  landmarks.insert(landmarks.begin(), landmark);


}

void mitk::NavigationDataXDofTo6DofFilter::PassThrough(unsigned int inputID, unsigned int outputID)
{

  // get the needed variables (input and output)
  const mitk::NavigationData *nd = this->GetInput(inputID);
  mitk::NavigationData *output = this->GetOutput(outputID);

  if (!nd || !output)
  {
      MITK_ERROR("NavigationDataXDofTo6DofFilter") << "Input and output must not be null.";
      return;
  }

  output->Graft(nd); // copy all information from input to output
  output->SetDataValid(nd->IsDataValid());

}

void mitk::NavigationDataXDofTo6DofFilter::UpdateListOfAllOutputs(){
  for (Landmark l : landmarks){
    if(!IsOutputInList(l.outputID)){m_OutputList.push_back(l.outputID);}
  }
}

bool mitk::NavigationDataXDofTo6DofFilter::IsOutputInList(int outputID){
  for(int id : m_OutputList){
    if(id==outputID){return true;}
  }
  return false;
}
