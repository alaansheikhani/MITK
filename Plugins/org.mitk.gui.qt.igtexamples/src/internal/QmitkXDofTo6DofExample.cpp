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



const std::string QmitkXDofTo6DofExample::VIEW_ID = "org.mitk.views.QmitkXDofTo6DofExample";

void QmitkXDofTo6DofExample::SetFocus()
{
}

QmitkXDofTo6DofExample::~QmitkXDofTo6DofExample()
{
}

void QmitkXDofTo6DofExample::CreateQtPartControl( QWidget *parent )
{

  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi( parent );

  // connect the widget items with the methods
  connect( m_Controls.butStart, SIGNAL(clicked()),
           this, SLOT(Start()) );
}


void QmitkXDofTo6DofExample::Start()
{
  MITK_INFO << "Test";
}
