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

#include "mitkTrackingDeviceSource.h"

#include "mitkTrackingDevice.h"
#include "mitkTrackingTool.h"

#include "mitkIGTTimeStamp.h"
#include "mitkIGTException.h"
#include "mitkIGTHardwareException.h"
#include <mitkNavigationDataLandmarkTransformFilter.h>

#include <vtkLandmarkTransform.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkPoints.h>

mitk::TrackingDeviceSource::TrackingDeviceSource()
  : mitk::NavigationDataSource(), m_TrackingDevice(nullptr)
{
}

mitk::TrackingDeviceSource::~TrackingDeviceSource()
{
  if (m_TrackingDevice.IsNotNull())
  {
    if (m_TrackingDevice->GetState() == mitk::TrackingDevice::Tracking)
    {
      this->StopTracking();
    }
    if (m_TrackingDevice->GetState() == mitk::TrackingDevice::Ready)
    {
      this->Disconnect();
    }
    m_TrackingDevice = nullptr;
  }
}

void mitk::TrackingDeviceSource::GenerateData()
{
  if (m_IsFrozen) {return;} //no update at all if device is frozen
  else if (m_TrackingDevice.IsNull()) {return;}

  if (m_TrackingDevice->GetToolCount() < 1)
    return;

  if (this->GetNumberOfIndexedOutputs() != m_TrackingDevice->GetToolCount()) // mismatch between tools and outputs. What should we do? Were tools added to the tracking device after SetTrackingDevice() was called?
  {
    //check this: TODO:
    ////this might happen if a tool is plugged into an aurora during tracking.
    //this->CreateOutputs();

    std::stringstream ss;
    ss << "mitk::TrackingDeviceSource: not enough outputs available for all tools. "
      << this->GetNumberOfOutputs() << " outputs available, but "
      << m_TrackingDevice->GetToolCount() << " tools available in the tracking device.";
    throw std::out_of_range(ss.str());
  }

  // Declare needed Variables for Landmark-Transformation

  //mitk::NavigationDataLandmarkTransformFilter::Pointer m_LandmarkTransformFilter = mitk::NavigationDataLandmarkTransformFilter::New();
  //mitk::PointSet::Pointer m_SourcePointSet = mitk::PointSet::New();
  //mitk::PointSet::Pointer m_TargetPointSet = mitk::PointSet::New();

  //mitk::PointSet::PointType pointonVirtual;
  //m_SourcePointSet->GetPointIfExists(0, &pointonVirtual);
  //mitk::NavigationData::Pointer m_virtualModelInputND = mitk::NavigationData::New();
  //m_virtualModelInputND->SetPosition(pointonVirtual);
  //m_virtualModelInputND->Modified();

  vtkSmartPointer<vtkPoints> sourcePoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> targetPoints = vtkSmartPointer<vtkPoints>::New();


  //Define and Insert SourcePoints
  // mitk::NavigationData::PositionType srcPt1;
  // mitk::Point3D srcPt1;
  // srcPt1[0] = 0;
  // srcPt1[1] = 20;
  // srcPt1[2] = 0;
  double srcPt1[3] = {0, 20, 0};
  double srcPt2[3] = {0, 0, 0};
  double srcPt3[3] = {20, 0, 0};
  double srcPt4[3] = {40, 0, 0};

  //m_SourcePointSet->InsertPoint(srcPt1);
  //m_SourcePointSet->InsertPoint(srcPt2);
  //m_SourcePointSet->InsertPoint(srcPt3);
  //m_SourcePointSet->InsertPoint(srcPt4);
  sourcePoints->InsertNextPoint(srcPt1);
  sourcePoints->InsertNextPoint(srcPt2);
  sourcePoints->InsertNextPoint(srcPt3);
  sourcePoints->InsertNextPoint(srcPt4);


  /* update outputs with tracking data from tools */
  unsigned int toolCount = m_TrackingDevice->GetToolCount();
  for (unsigned int i = 0; i < toolCount; ++i)
  {
    mitk::NavigationData* nd = this->GetOutput(i);
    assert(nd);
    mitk::TrackingTool* t = m_TrackingDevice->GetTool(i);
    assert(t);

    if ((t->IsEnabled() == false) || (t->IsDataValid() == false))
    {
      nd->SetDataValid(false);
      continue;
    }
    nd->SetDataValid(true);
    mitk::NavigationData::PositionType p;
    mitk::NavigationData::OrientationType o;
    t->GetPosition(p);
    nd->SetPosition(p);
    t->GetOrientation(o);
    nd->SetOrientation(o);
    nd->SetOrientationAccuracy(t->GetTrackingError());
    nd->SetPositionAccuracy(t->GetTrackingError());
    nd->SetIGTTimeStamp(t->GetIGTTimeStamp());

    if (i == 0)
    {
      double targetPoint1[3] = {0, 0, 0};
      mitk::NavigationData::Pointer trgPt1 = mitk::NavigationData::New();
      trgPt1->SetPosition(targetPoint1);
      trgPt1->Compose(nd);

      double targetPoint2[3] = {0, 0, 20};
      mitk::NavigationData::Pointer trgPt2 = mitk::NavigationData::New();
      trgPt2->SetPosition(targetPoint2);
      trgPt2->Compose(nd);

      //m_TargetPointSet->InsertPoint(trgPt1->GetPosition());
      //m_TargetPointSet->InsertPoint(trgPt2->GetPosition());
      mitk::Point3D t1 = trgPt1->GetPosition();
      double targetPointComp[3] = {t1[0], t1[1], t1[2]};
      targetPoints->InsertNextPoint(targetPointComp);

      mitk::Point3D t2 = trgPt2->GetPosition();
      double targetPointComp2[3] = {t2[0], t2[1], t2[2]};
      targetPoints->InsertNextPoint(targetPointComp2);
    }
    if (i == 1)
    {
      double targetPoint3[3] = {0, 0, 0};
      mitk::NavigationData::Pointer trgPt3 = mitk::NavigationData::New();
      trgPt3->SetPosition(targetPoint3);
      trgPt3->Compose(nd);
      //trgPt1->Compose(this->GetOutput(0));


      mitk::Point3D t3 = trgPt3->GetPosition();
      double targetPointComp3[3] = {t3[0], t3[1], t3[2]};
      targetPoints->InsertNextPoint(targetPointComp3);

      double targetPoint4[3] = {0, 0, 20};

      mitk::NavigationData::Pointer trgPt4 = mitk::NavigationData::New();
      trgPt4->SetPosition(targetPoint4);
      trgPt4->Compose(nd);

      //m_TargetPointSet->InsertPoint(trgPt3->GetPosition());
      //m_TargetPointSet->InsertPoint(trgPt4->GetPosition());

      mitk::Point3D t4 = trgPt4->GetPosition();
      double targetPointComp4[4] = {t4[0], t4[1], t4[2]};
      targetPoints->InsertNextPoint(targetPointComp4);

      //m_LandmarkTransformFilter->SetSourceLandmarks(m_SourcePointSet);
      //m_LandmarkTransformFilter->SetTargetLandmarks(m_TargetPointSet);
      //m_LandmarkTransformFilter->SetInput(0, m_virtualModelInputND);
      //m_LandmarkTransformFilter->Update();
      //m_LandmarkTransformFilter->GetOutput(0);

      // compute transform
      vtkSmartPointer<vtkLandmarkTransform> transform = vtkSmartPointer<vtkLandmarkTransform>::New();
      transform->SetSourceLandmarks(sourcePoints);
      transform->SetTargetLandmarks(targetPoints);
      transform->SetModeToRigidBody();
      transform->Modified();
      transform->Update();

      // convert from vtk to itk data types
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
      // create affine transform 3D surface
      mitk::AffineTransform3D::Pointer mitkTransform = mitk::AffineTransform3D::New();
      mitkTransform->SetMatrix(rotationDouble);
      mitkTransform->SetOffset(translationDouble);
      //#############################################################################################

      //############### object is transformed ##########################
      // save transform
      mitk::NavigationData::Pointer transformedND = mitk::NavigationData::New(mitkTransform); // this is stored in a member because it is needed for permanent registration later on
      transformedND->SetName("Test6D");
      transformedND->SetDataValid(true);
      this->GetOutput(i)->Graft(transformedND);

     //mitk::NavigationData::Pointer transformedND = mitk::NavigationData::New();
     //nd = transformedND;
     //nd->SetPosition(transformedND->GetPosition());
     //nd->SetOrientation(transformedND->GetOrientation());

     }

    //nd->SetPosition(p);

    //mitk::NavigationData::OrientationType o;
    //t->GetOrientation(o);
    //nd->SetOrientation(o);
    ////nd->SetOrientationAccuracy(t->GetTrackingError());
    ////nd->SetPositionAccuracy(t->GetTrackingError());
    ////nd->SetIGTTimeStamp(t->GetIGTTimeStamp());

    //for backward compatibility: check if the timestamp was set, if not create a default timestamp
    //if (nd->GetIGTTimeStamp()==0) nd->SetIGTTimeStamp(mitk::IGTTimeStamp::GetInstance()->GetElapsed());

  }

  //mitk::NavigationData::Pointer test1 = mitk::NavigationData::New();
  //mitk::NavigationData::Pointer test2 = mitk::NavigationData::New();
  //test1->SetName("Test");
  //test1->SetDataValid(true);
  //this->GetOutput(0)->Graft(test1);
  //this->GetOutput(1)->Graft(test2);
  //double point[3] = {0, 0, 0};
  /*
  test1->SetPosition(point);
  test1->SetDataValid(true);
  test2->SetPosition(point);

  this->SetNthOutput(0, test1);
  this->SetNthOutput(1, test2);*/
}


void mitk::TrackingDeviceSource::SetTrackingDevice( mitk::TrackingDevice* td )
{
  MITK_DEBUG << "Setting TrackingDevice to " << td;
  if (this->m_TrackingDevice.GetPointer() != td)
  {
    this->m_TrackingDevice = td;
    this->CreateOutputs();
    std::stringstream name; // create a human readable name for the source
    name << td->GetData().Model << " Tracking Source";
    this->SetName(name.str());
  }
}

void mitk::TrackingDeviceSource::CreateOutputs(){
  //if outputs are set then delete them
  if (this->GetNumberOfOutputs() > 0)
  {
    for (int numOP = this->GetNumberOfOutputs() -1; numOP >= 0; numOP--)
      this->RemoveOutput(numOP);
    this->Modified();
  }

  //fill the outputs if a valid tracking device is set
  if (m_TrackingDevice.IsNull())
    return;

  this->SetNumberOfIndexedOutputs(m_TrackingDevice->GetToolCount());  // create outputs for all tools
  unsigned int numberOfOutputs = this->GetNumberOfIndexedOutputs();
  MITK_DEBUG << "Number of tools at start of method CreateOutputs(): " << m_TrackingDevice->GetToolCount();
  MITK_DEBUG << "Number of outputs at start of method CreateOutputs(): " << numberOfOutputs;
  for (unsigned int idx = 0; idx < m_TrackingDevice->GetToolCount(); ++idx)
  {
    if (this->GetOutput(idx) == nullptr)
    {
      DataObjectPointer newOutput = this->MakeOutput(idx);
      static_cast<mitk::NavigationData*>(newOutput.GetPointer())->SetName(m_TrackingDevice->GetTool(idx)->GetToolName()); // set NavigationData name to ToolName
      this->SetNthOutput(idx, newOutput);
      this->Modified();
    }
  }
}

void mitk::TrackingDeviceSource::Connect()
{
  if (m_TrackingDevice.IsNull())
    throw std::invalid_argument("mitk::TrackingDeviceSource: No tracking device set");
  if (this->IsConnected())
    return;
  try
  {
    //Try to open the connection. If it didn't work (fals is returned from OpenConnection by the tracking device), throw an exception.
    if (!m_TrackingDevice->OpenConnection())
    {
      mitkThrowException(mitk::IGTHardwareException) << "Could not open connection.";
    }
  }
  catch (mitk::IGTException &e)
  {
    throw std::runtime_error(std::string("mitk::TrackingDeviceSource: Could not open connection to tracking device. Error: ") + e.GetDescription());
  }
}

void mitk::TrackingDeviceSource::StartTracking()
{
  if (m_TrackingDevice.IsNull())
    throw std::invalid_argument("mitk::TrackingDeviceSource: No tracking device set");
  if (m_TrackingDevice->GetState() == mitk::TrackingDevice::Tracking)
    return;
  if (m_TrackingDevice->StartTracking() == false)
    throw std::runtime_error("mitk::TrackingDeviceSource: Could not start tracking");
}

void mitk::TrackingDeviceSource::Disconnect()
{
  if (m_TrackingDevice.IsNull())
    throw std::invalid_argument("mitk::TrackingDeviceSource: No tracking device set");
  if (m_TrackingDevice->CloseConnection() == false)
    throw std::runtime_error("mitk::TrackingDeviceSource: Could not close connection to tracking device");
}

void mitk::TrackingDeviceSource::StopTracking()
{
  if (m_TrackingDevice.IsNull())
    throw std::invalid_argument("mitk::TrackingDeviceSource: No tracking device set");
  if (m_TrackingDevice->StopTracking() == false)
    throw std::runtime_error("mitk::TrackingDeviceSource: Could not stop tracking");
}

void mitk::TrackingDeviceSource::UpdateOutputInformation()
{
  if(this->GetTrackingDevice()->GetToolCount() != this->GetNumberOfIndexedOutputs())
    this->CreateOutputs();

  this->Modified();  // make sure that we need to be updated
  Superclass::UpdateOutputInformation();
}

//unsigned int mitk::TrackingDeviceSource::GetToolCount()
//{
//  if (m_TrackingDevice)
//    return m_TrackingDevice->GetToolCount();
//  return 0;
//}

bool mitk::TrackingDeviceSource::IsConnected()
{
  if (m_TrackingDevice.IsNull())
    return false;

  return (m_TrackingDevice->GetState() == mitk::TrackingDevice::Ready) || (m_TrackingDevice->GetState() == mitk::TrackingDevice::Tracking);
}

bool mitk::TrackingDeviceSource::IsTracking()
{
  if (m_TrackingDevice.IsNull())
    return false;

  return m_TrackingDevice->GetState() == mitk::TrackingDevice::Tracking;
}
