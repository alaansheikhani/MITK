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

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkRenderWindow.h"

#include "QmitkXDofTo6DofExample.h"
#include "mitkNavigationDataObjectVisualizationFilter.h"
#include <vtkMarchingCubes.h>
#include <mitkSurface.h>
#include <mitkCone.h>

const std::string QmitkXDofTo6DofExample::VIEW_ID = "org.mitk.views.QmitkXDofTo6DofExample";

void QmitkXDofTo6DofExample::SetFocus() {}

QmitkXDofTo6DofExample::~QmitkXDofTo6DofExample() {}

void QmitkXDofTo6DofExample::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  m_Timer = new QTimer(this);

  // connect the widget items with the methods
  connect(m_Controls.butStart, SIGNAL(clicked()), this, SLOT(AddLandmark()));
  connect(m_Timer, SIGNAL(timeout()), this, SLOT(Update()));
  // connect(m_Controls.but)
  // m_TrackingDeviceSource->SetTrackingDevice()
}

void QmitkXDofTo6DofExample::Update()
{
    //hier Vis Filter updaten
  MITK_INFO << "Update!";
}

void QmitkXDofTo6DofExample::AddLandmark()
{
  static bool isFirstTime = true;
  if (isFirstTime)
  {
    // this->CreatePipeline();
    isFirstTime = false;
  }
  /*
    m_Timer.setInterval(this->m_Controls.visualizationUpdateRateSpinBox->value());
    m_Timer.start();
    m_Controls.*/

  MITK_INFO << "Test";
  m_Source = m_Controls.NavigationDataSourceSelectionWidget->GetSelectedNavigationDataSource();
  // m_Controls.NavigationDataSourceSelectionWidget->GetSelectedNavigationDataSource()
  int inputID = m_Controls.NavigationDataSourceSelectionWidget->GetSelectedToolID();

  double source1 = m_Controls.source1->value();
  double source2 = m_Controls.source2->value();
  double source3 = m_Controls.source3->value();

  double target1 = m_Controls.target1->value();
  double target2 = m_Controls.target2->value();
  double target3 = m_Controls.target3->value();

  int outputID = m_Controls.outputID->value();



  //Surface nötig für Data Node

  mitk::Cone::Pointer cone = mitk::Cone::New();
  //mitk::Surface::Pointer surface = mitk::Surface::New();


  mitk::DataNode::Pointer representationObject = mitk::DataNode::New();
  representationObject->SetName("X Dof to 6 Dof Example Object");
  //set name, set data (z.B. Cone)
  representationObject->SetData(cone);
  this->GetDataStorage()->Add(representationObject);
  // 1. Parameter 0, wenn nur ein Tool
  m_VisFilter->SetRepresentationObject(0, representationObject->GetData());
}

void QmitkXDofTo6DofExample::StartUpdating() {

    //Hier Pipeline aufbauen

    m_XDofTo6DofFilter = mitk::NavigationDataXDofTo6DofFilter::New();
    m_XDofTo6DofFilter->ConnectTo(m_Source);
    m_VisFilter = mitk::NavigationDataObjectVisualizationFilter::New();
    m_VisFilter->ConnectTo(m_XDofTo6DofFilter);

    m_Timer->start(200);

}

// void QmitkXDofTo6DofExample::CreatePipeline()
//{
//  // create a visualization filter
//  m_VisFilter = mitk::NavigationDataObjectVisualizationFilter::New();
//
//}
