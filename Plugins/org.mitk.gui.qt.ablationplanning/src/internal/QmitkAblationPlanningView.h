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


#ifndef QMITKABLATIONPLANNINGVIEW_H
#define QMITKABLATIONPLANNINGVIEW_H

#include <berryISelectionListener.h>

#include <QmitkAbstractView.h>

#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateOr.h>

#include "ui_QmitkAblationPlanningViewControls.h"




/**
  \brief QmitkAblationPlanningView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class QmitkAblationPlanningView : public QmitkAbstractView
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:

  // a type for handling lists of DataNodes
  typedef std::vector<mitk::DataNode*> NodeList;

  typedef std::map<mitk::DataNode*, unsigned long> NodeTagMapType;


  QmitkAblationPlanningView();
  virtual ~QmitkAblationPlanningView();

  /*!
  \brief Invoked when the DataManager selection changed
  */
  virtual void OnSelectionChanged(mitk::DataNode* node);
  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part,
    const QList<mitk::DataNode::Pointer>& nodes) override;

  static const std::string VIEW_ID;

protected:
  virtual void CreateQtPartControl(QWidget *parent) override;
  virtual void SetFocus() override;

  //void ResetMouseCursor();
  //void SetMouseCursor(const us::ModuleResource&, int hotspotX, int hotspotY);

  void NodeRemoved(const mitk::DataNode* node) override;

  void NodeAdded(const mitk::DataNode *node) override;

  bool CheckForSameGeometry(const mitk::DataNode*, const mitk::DataNode*) const;

protected slots:
  void OnSegmentationComboBoxSelectionChanged(const mitk::DataNode* node);
  void OnVisiblePropertyChanged();
  void OnBinaryPropertyChanged();
  void OnAblationStatingPointPushButtonClicked();


private:
  Ui::QmitkAblationPlanningViewControls m_Controls;

  mitk::NodePredicateAnd::Pointer m_IsOfTypeImagePredicate;
  mitk::NodePredicateOr::Pointer m_IsASegmentationImagePredicate;
  mitk::NodePredicateAnd::Pointer m_IsAPatientImagePredicate;

  NodeTagMapType  m_WorkingDataObserverTags;
  NodeTagMapType  m_BinaryPropertyObserverTags;

  unsigned long m_VisibilityChangedObserverTag;
  bool m_MouseCursorSet;
  bool m_DataSelectionChanged;

  mitk::Point3D m_AblationStartingPositionInWorldCoordinates;
  itk::Index<3> m_AblationStartingPositionIndexCoordinates;
  bool m_AblationStartingPositionValid;
};

#endif // QMITKABLATIONPLANNINGVIEW_H
