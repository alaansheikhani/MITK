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

#ifndef QmitkXDofTo6DofExample_h
#define QmitkXDofTo6DofExample_h

#include <berryISelectionListener.h>

#include <QmitkAbstractView.h>

#include "ui_QmitkXDofTo6DofExample.h"

#include "mitkNavigationDataObjectVisualizationFilter.h"
#include "mitkTrackingDeviceSource.h"
#include "mitkNavigationDataXDofTo6DofFilter.h"
#include <QTimer>

/**
  \brief OpenIGTLinkExample

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class QmitkXDofTo6DofExample : public QmitkAbstractView
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:
  ~QmitkXDofTo6DofExample();

  static const std::string VIEW_ID;

protected slots:

  void AddLandmark();
  void StartUpdating();
  void Update();
  //void CreatePipeline();

protected:
  Ui::QmitkXDofTo6DofExampleControls m_Controls;

  virtual void CreateQtPartControl(QWidget *parent) override;

  virtual void SetFocus() override;

  mitk::NavigationDataXDofTo6DofFilter::Pointer m_XDofTo6DofFilter;
  mitk::NavigationDataObjectVisualizationFilter::Pointer m_VisFilter;
  mitk::NavigationDataSource::Pointer m_Source;
  QTimer* m_Timer;
};

#endif // QmitkXDofTo6DofExample_h
