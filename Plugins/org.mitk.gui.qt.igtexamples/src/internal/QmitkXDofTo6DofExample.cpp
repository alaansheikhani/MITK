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
  connect(m_Controls.butStartUpdating, SIGNAL(clicked()), this, SLOT(StartUpdating()));
}

void QmitkXDofTo6DofExample::Update()
{
  m_VisFilter->Update();
}

void QmitkXDofTo6DofExample::AddLandmark()
{
  MITK_INFO << "Add Landmark";
  unsigned int inputID = m_Controls.NavigationDataSourceSelectionWidget->GetSelectedToolID();

  double source1 = m_Controls.source1->value();
  double source2 = m_Controls.source2->value();
  double source3 = m_Controls.source3->value();
  mitk::Point3D source;
  mitk::FillVector3D(source,source1,source2,source3);

  double target1 = m_Controls.target1->value();
  double target2 = m_Controls.target2->value();
  double target3 = m_Controls.target3->value();
  mitk::Point3D target;
  mitk::FillVector3D(target,target1,target2,target3);

  unsigned int outputID = m_Controls.outputID->value();

  m_XDofTo6DofFilter->AddLandmarkFor6DoF(source,target,inputID,outputID);
}

void QmitkXDofTo6DofExample::StartUpdating() {
    m_Source = m_Controls.NavigationDataSourceSelectionWidget->GetSelectedNavigationDataSource();
    if(m_Source.IsNull()){
      MITK_WARN << "Cannot start updating, no source selected!";
      return;
    }

    //Connect / build up pipeline
    m_XDofTo6DofFilter = mitk::NavigationDataXDofTo6DofFilter::New();
    m_XDofTo6DofFilter->ConnectTo(m_Source);

    m_VisFilter = mitk::NavigationDataObjectVisualizationFilter::New();
    m_VisFilter->ConnectTo(m_XDofTo6DofFilter);
    for(int id=0;id<m_VisFilter->GetNumberOfInputs();id++){
       MITK_INFO << "Add representation object for input " << id;
       mitk::Cone::Pointer cone = mitk::Cone::New();
       mitk::DataNode::Pointer representationObject = mitk::DataNode::New();
       representationObject->SetData(cone);
       representationObject->SetName("X Dof to 6 Dof Example Object");
       this->GetDataStorage()->Add(representationObject);
       m_VisFilter->SetRepresentationObject(id, representationObject->GetData());
    }

    //start updating pipeline
    m_Timer->start(200);
    MITK_INFO << "Started updating X DoF to 6 DoF example IGT pipeline!";
}
