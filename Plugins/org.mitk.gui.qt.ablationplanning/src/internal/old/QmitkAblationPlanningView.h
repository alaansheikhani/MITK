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

#include <mitkImage.h>
#include "mitkAblationPlan.h"
#include "mitkAblationPlanningAlgorithm.h"
#include "mitkAblationUtils.h"
#include "mitkAblationZone.h"

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
  typedef std::vector<mitk::DataNode *> NodeList;

  typedef std::map<mitk::DataNode *, unsigned long> NodeTagMapType;

  QmitkAblationPlanningView();
  virtual ~QmitkAblationPlanningView();

  /*!
  \brief Invoked when the DataManager selection changed
  */
  virtual void OnSelectionChanged(mitk::DataNode *node);
  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part,
                                  const QList<mitk::DataNode::Pointer> &nodes) override;

  void UnsetSegmentationImageGeometry();
  void SetSegmentationImageGeometryInformation(mitk::Image *image);

  static const std::string VIEW_ID;

protected:
  virtual void CreateQtPartControl(QWidget *parent) override;
  virtual void SetFocus() override;

  void NodeRemoved(const mitk::DataNode *node) override;

  void NodeAdded(const mitk::DataNode *node) override;

  bool CheckForSameGeometry(const mitk::DataNode *, const mitk::DataNode *) const;

  void CopyTemporaryAblationZoneDistribution();

  void CreateSpheresOfAblationVolumes();

  void DeleteAllSpheres();

  void CalculateAblationStatistics();

protected slots:
  void OnSegmentationComboBoxSelectionChanged(const mitk::DataNode *node);
  void OnVisiblePropertyChanged();
  void OnBinaryPropertyChanged();
  void OnCalculateSafetyMargin();
  void OnAblationStartingPointPushButtonClicked();
  void OnCalculateAblationZonesPushButtonClicked();

private:
  Ui::QmitkAblationPlanningViewControls m_Controls;

  mitk::NodePredicateAnd::Pointer m_IsOfTypeImagePredicate;
  mitk::NodePredicateOr::Pointer m_IsASegmentationImagePredicate;
  mitk::NodePredicateAnd::Pointer m_IsAPatientImagePredicate;

  NodeTagMapType m_WorkingDataObserverTags;
  NodeTagMapType m_BinaryPropertyObserverTags;

  unsigned long m_VisibilityChangedObserverTag;
  bool m_MouseCursorSet;
  bool m_DataSelectionChanged;

  mitk::Point3D m_AblationStartingPositionInWorldCoordinates; //entf
  itk::Index<3> m_AblationStartingPositionIndexCoordinates; //entf
  mitk::Point3D m_TempAblationStartingPositionInWorldCoordinates; //entf
  itk::Index<3> m_TempAblationStartingPositionIndexCoordinates; //entf

  bool m_ManualAblationStartingPositionSet; //entf
  double m_MaxAblationRadius; // Maximal ablation radius //entf
  double m_AblationRadius; // Desired ablation radius //entf
  double m_MinAblationRadius; // Minimal ablation radius //entf
  mitk::Image::Pointer m_SegmentationImage;
  mitk::AblationPlan::Pointer m_AblationPlan;
  mitk::AblationPlanningAlgorithm::Pointer m_PlanningAlgo;
  /*!
  * \brief Final vector storing the index coordinates of all circle centers of the ablation zones
    after the calculation of the best ablation zone distribution.
  */
  std::vector<mitk::AblationZone> m_AblationZones;

  /*!
   * \brief Temporary vector storing the index coordinates of all circle centers of the ablation zones
   * when calculating the best ablation zone distribution.
   */
  std::vector<mitk::AblationZone> m_TempAblationZones; //entf?

  /*!
   * \brief Vector storing the index coordinates of all circle centers of the ablation zones,
   * which are finally processed after the calculation of the best ablation zone distribution.
   * This means: All 12 direct neighbour ablation zones are checked for remaining non-ablated tumor issue.
   */
  std::vector<mitk::AblationZone> m_AblationZonesProcessed;

  /*!
   * \brief Temporary vector storing the index coordinates of all circle centers of the ablation zones,
   * which are finally processed. This means: All 12 direct neighbour ablation zones are checked
   * for remaining non-ablated tumor issue.
   */
  std::vector<mitk::AblationZone> m_TempAblationZonesProcessed;

  /*!
   * \brief Vector storing the index coordinates of all pixels, which are tumor tissue or
   * safety margin.
   */
  std::vector<itk::Index<3>> m_TumorTissueSafetyMarginIndices;

  /*!
  \brief The 3D dimension of the segmentation image given in index size.
  */
  mitk::Vector3D m_ImageDimension;
  mitk::Vector3D m_ImageSpacing;

  /** Holds a point set with the centers of all ablation zones */
  mitk::DataNode::Pointer m_AblationCentersNode;

  bool m_AblationCalculationMade;
};

#endif // QMITKABLATIONPLANNINGVIEW_H
