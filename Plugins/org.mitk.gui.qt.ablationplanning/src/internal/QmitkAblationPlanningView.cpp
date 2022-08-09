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
#include "QmitkAblationPlanningView.h"

// Qt
#include <QMessageBox>
#include <qfiledialog.h>

// mitk
#include "mitkProperties.h"
#include <ctime>
#include <mitkImage.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkPointSet.h>
#include <mitkSurface.h>
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>

const std::string QmitkAblationPlanningView::VIEW_ID = "org.mitk.views.ablationplanning";
const static short ABLATION_VALUE = 2;
const static short TUMOR_NOT_YET_ABLATED = 1;
const static short NO_TUMOR_ISSUE = 0;
const static short SAFETY_MARGIN = 256;
const static unsigned short BIT_OPERATION_ELIMINATE_TUMOR_SAFETY_MARGIN = 65278; // = 11111110 11111110

//=====================Konstruktor/Destruktor===================================
QmitkAblationPlanningView::QmitkAblationPlanningView()
  : m_MouseCursorSet(false),
    m_DataSelectionChanged(false),
    m_AblationStartingPositionInWorldCoordinates(),
    m_AblationStartingPositionIndexCoordinates(),
    m_ManualAblationStartingPositionSet(false),
    m_AblationCalculationMade(false),
    m_AblationCentersNode(mitk::DataNode::New()),
    m_PlanningAlgo(mitk::AblationPlanningAlgorithm::New()),
    m_PlanLogger(mitk::AblationPlanningLogging::New())
{
  this->UnsetSegmentationImageGeometry();

  mitk::TNodePredicateDataType<mitk::Image>::Pointer isImage = mitk::TNodePredicateDataType<mitk::Image>::New();
  auto isSegmentation = mitk::NodePredicateDataType::New("Segment");

  mitk::NodePredicateOr::Pointer validImages = mitk::NodePredicateOr::New();
  validImages->AddPredicate(mitk::NodePredicateAnd::New(isImage, mitk::NodePredicateNot::New(isSegmentation)));

  mitk::NodePredicateNot::Pointer isNotAHelperObject =
    mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true)));

  m_IsOfTypeImagePredicate = mitk::NodePredicateAnd::New(validImages, isNotAHelperObject);

  mitk::NodePredicateProperty::Pointer isBinaryPredicate =
    mitk::NodePredicateProperty::New("binary", mitk::BoolProperty::New(true));
  mitk::NodePredicateNot::Pointer isNotBinaryPredicate = mitk::NodePredicateNot::New(isBinaryPredicate);

  mitk::NodePredicateAnd::Pointer isABinaryImagePredicate =
    mitk::NodePredicateAnd::New(m_IsOfTypeImagePredicate, isBinaryPredicate);
  mitk::NodePredicateAnd::Pointer isNotABinaryImagePredicate =
    mitk::NodePredicateAnd::New(m_IsOfTypeImagePredicate, isNotBinaryPredicate);

  m_IsASegmentationImagePredicate =
    mitk::NodePredicateOr::New(isABinaryImagePredicate, mitk::TNodePredicateDataType<mitk::LabelSetImage>::New());
  m_IsAPatientImagePredicate = mitk::NodePredicateAnd::New(
    isNotABinaryImagePredicate, mitk::NodePredicateNot::New(mitk::TNodePredicateDataType<mitk::LabelSetImage>::New()));

  // to initialize with the first image, if nothing was selected yet...
  if (this->GetDataStorage()->GetNode(m_IsASegmentationImagePredicate) != nullptr)
  {
    m_SegmentationImage =
      dynamic_cast<mitk::Image *>(this->GetDataStorage()->GetNode(m_IsASegmentationImagePredicate)->GetData());
    SetSegmentationImageGeometryInformation(m_SegmentationImage);
    OnSelectionChanged(this->GetDataStorage()->GetNode(m_IsASegmentationImagePredicate));
  }
}

QmitkAblationPlanningView::~QmitkAblationPlanningView()
{
  // removing all observers
  for (NodeTagMapType::iterator dataIter = m_WorkingDataObserverTags.begin();
       dataIter != m_WorkingDataObserverTags.end();
       ++dataIter)
  {
    (*dataIter).first->GetProperty("visible")->RemoveObserver((*dataIter).second);
  }
  m_WorkingDataObserverTags.clear();

  for (NodeTagMapType::iterator dataIter = m_BinaryPropertyObserverTags.begin();
       dataIter != m_BinaryPropertyObserverTags.end();
       ++dataIter)
  {
    (*dataIter).first->GetProperty("binary")->RemoveObserver((*dataIter).second);
  }
  m_BinaryPropertyObserverTags.clear();
}

void QmitkAblationPlanningView::OnSelectionChanged(mitk::DataNode *node)
{
  MITK_DEBUG << "OnSelectionChanged()";
  berry::IWorkbenchPart::Pointer nullPart;
  QList<mitk::DataNode::Pointer> nodes;
  nodes.push_back(node);
  this->OnSelectionChanged(nullPart, nodes);
}

void QmitkAblationPlanningView::OnSelectionChanged(berry::IWorkbenchPart::Pointer part,
                                                   const QList<mitk::DataNode::Pointer> &nodes)
{
  MITK_DEBUG << "OnSelectionChanged()";

  if (nodes.size() == 1)
  {
    mitk::DataNode::Pointer selectedNode = nodes.at(0);
    if (selectedNode.IsNull())
    {
      return;
    }

    mitk::Image::Pointer selectedImage = dynamic_cast<mitk::Image *>(selectedNode->GetData());
    if (selectedImage.IsNull())
    {
      return;
    }

    if (m_IsASegmentationImagePredicate->CheckNode(selectedNode))
    {
      // if a segmentation is selected find a possible patient image
      mitk::DataStorage::SetOfObjects::ConstPointer sources =
        GetDataStorage()->GetSources(selectedNode, m_IsAPatientImagePredicate);
      mitk::DataNode::Pointer refNode;
      if (sources->Size() != 0)
      {
        // found one or more sources - use the first one
        refNode = sources->ElementAt(0);
        refNode->SetVisibility(true);
        selectedNode->SetVisibility(true);
      }
      // mitk::RenderingManager::GetInstance()->InitializeViews(selectedNode->GetData()->GetTimeGeometry(),
      // mitk::RenderingManager::REQUEST_UPDATE_ALL, true);
    }
    else
    {
      MITK_WARN << "SelectedNode is no segmentation node";
    }
  }
}

void QmitkAblationPlanningView::UnsetSegmentationImageGeometry()
{
  m_ImageDimension[0] = 0;
  m_ImageDimension[1] = 0;
  m_ImageDimension[2] = 0;

  m_ImageSpacing[0] = 1;
  m_ImageSpacing[1] = 1;
  m_ImageSpacing[2] = 1;
}

void QmitkAblationPlanningView::SetSegmentationImageGeometryInformation(mitk::Image *image)
{
  m_ImageDimension[0] = image->GetDimension(0);
  m_ImageDimension[1] = image->GetDimension(1);
  m_ImageDimension[2] = image->GetDimension(2);

  m_ImageSpacing[0] = image->GetGeometry()->GetSpacing()[0];
  m_ImageSpacing[1] = image->GetGeometry()->GetSpacing()[1];
  m_ImageSpacing[2] = image->GetGeometry()->GetSpacing()[2];
}

//==============================================================================

void QmitkAblationPlanningView::SetFocus()
{
  m_Controls.segmentationComboBox->setFocus();
}

void QmitkAblationPlanningView::NodeRemoved(const mitk::DataNode *node)
{
  MITK_DEBUG << "NodeRemoved()";
  if (m_IsASegmentationImagePredicate->CheckNode(node))
  {
    // First of all remove all possible contour markers of the segmentation
    mitk::DataStorage::SetOfObjects::ConstPointer allContourMarkers = this->GetDataStorage()->GetDerivations(
      node, mitk::NodePredicateProperty::New("isContourMarker", mitk::BoolProperty::New(true)));

    for (mitk::DataStorage::SetOfObjects::ConstIterator it = allContourMarkers->Begin(); it != allContourMarkers->End();
         ++it)
    {
      std::string nodeName = node->GetName();
      unsigned int t = nodeName.find_last_of(" ");
      unsigned int id = atof(nodeName.substr(t + 1).c_str()) - 1;

      this->GetDataStorage()->Remove(it->Value());
    }

    mitk::Image *image = dynamic_cast<mitk::Image *>(node->GetData());
  }
  mitk::DataNode *tempNode = const_cast<mitk::DataNode *>(node);
  // Since the binary property could be changed during runtime by the user
  if (m_IsOfTypeImagePredicate->CheckNode(node))
  {
    node->GetProperty("visible")->RemoveObserver(m_WorkingDataObserverTags[tempNode]);
    m_WorkingDataObserverTags.erase(tempNode);
    node->GetProperty("binary")->RemoveObserver(m_BinaryPropertyObserverTags[tempNode]);
    m_BinaryPropertyObserverTags.erase(tempNode);
  }
}

void QmitkAblationPlanningView::NodeAdded(const mitk::DataNode *node)
{
  MITK_DEBUG << "NodeAdded()";
  if (!m_IsOfTypeImagePredicate->CheckNode(node))
  {
    return;
  }

  itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command =
    itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
  command->SetCallbackFunction(this, &QmitkAblationPlanningView::OnVisiblePropertyChanged);
  m_WorkingDataObserverTags.insert(std::pair<mitk::DataNode *, unsigned long>(
    const_cast<mitk::DataNode *>(node), node->GetProperty("visible")->AddObserver(itk::ModifiedEvent(), command)));

  itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command2 =
    itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
  command2->SetCallbackFunction(this, &QmitkAblationPlanningView::OnBinaryPropertyChanged);
  m_BinaryPropertyObserverTags.insert(std::pair<mitk::DataNode *, unsigned long>(
    const_cast<mitk::DataNode *>(node), node->GetProperty("binary")->AddObserver(itk::ModifiedEvent(), command2)));
}

bool QmitkAblationPlanningView::CheckForSameGeometry(const mitk::DataNode *node1, const mitk::DataNode *node2) const
{
  bool isSameGeometry(true);

  mitk::Image *image1 = dynamic_cast<mitk::Image *>(node1->GetData());
  mitk::Image *image2 = dynamic_cast<mitk::Image *>(node2->GetData());
  if (image1 && image2)
  {
    mitk::BaseGeometry *geo1 = image1->GetGeometry();
    mitk::BaseGeometry *geo2 = image2->GetGeometry();

    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetOrigin(), geo2->GetOrigin());
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(0), geo2->GetExtent(0));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(1), geo2->GetExtent(1));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(2), geo2->GetExtent(2));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetSpacing(), geo2->GetSpacing());
    isSameGeometry = isSameGeometry && mitk::MatrixEqualElementWise(geo1->GetIndexToWorldTransform()->GetMatrix(),
                                                                    geo2->GetIndexToWorldTransform()->GetMatrix());

    return isSameGeometry;
  }
  else
  {
    return false;
  }
}

void QmitkAblationPlanningView::CopyTemporaryAblationZoneDistribution()
{
  if (m_AblationZonesProcessed.size() == 0 || m_AblationZonesProcessed.size() > m_TempAblationZonesProcessed.size())
  {
    MITK_INFO << "Reduced the number of ablation zones from: " << m_AblationZonesProcessed.size()
              << " to: " << m_TempAblationZonesProcessed.size() << "========================================";

    m_AblationZones.clear();
    m_AblationZonesProcessed.clear();
    m_AblationZones = m_TempAblationZones;
    m_AblationZonesProcessed = m_TempAblationZonesProcessed;
    m_TempAblationZones.clear();
    m_TempAblationZonesProcessed.clear();
    m_AblationStartingPositionIndexCoordinates = m_TempAblationStartingPositionIndexCoordinates;
    m_AblationStartingPositionInWorldCoordinates = m_TempAblationStartingPositionInWorldCoordinates;
  }
  else
  {
    m_TempAblationZones.clear();
    m_TempAblationZonesProcessed.clear();
  }
}

void QmitkAblationPlanningView::CreateSpheresOfAblationVolumes()
{
  mitk::PointSet::Pointer centerPoints = mitk::PointSet::New();
  // get tumor COG
  mitk::PointSet::Pointer COG = AblationUtils::CalculateCOGTargetPoints(selectedSurface);
  std::vector<mitk::PointSet::Pointer> zoneCenters;
  for (int index = 0; index < m_AblationPlan->GetNumberOfZones(); ++index)
  {
    mitk::DataNode::Pointer m_DataNode = mitk::DataNode::New();

    mitk::Surface::Pointer mySphere = mitk::Surface::New();

    vtkSmartPointer<vtkSphereSource> vtkSphere = vtkSmartPointer<vtkSphereSource>::New();
    vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    mitk::Point3D centerInWorldCoordinates;

    m_AblationPlan->GetSegmentationImage()->GetGeometry()->IndexToWorld(
      m_AblationPlan->GetAblationZone(index)->indexCenter, centerInWorldCoordinates);

    // Center
    vtkSphere->SetRadius(m_AblationPlan->GetAblationZone(index)->radius);
    vtkSphere->SetPhiResolution(40);
    vtkSphere->SetThetaResolution(40);
    vtkSphere->SetCenter(centerInWorldCoordinates[0], centerInWorldCoordinates[1], centerInWorldCoordinates[2]);
    vtkSphere->Update();

    appendPolyData->AddInputData(vtkSphere->GetOutput());

    mySphere->SetVtkPolyData(vtkSphere->GetOutput());

    // Add Node
    m_DataNode->SetData(mySphere);
    QString name = QString("Kugel_%1").arg(index + 1);
    m_DataNode->SetName(name.toStdString());
    m_DataNode->SetOpacity(0.3);
    this->GetDataStorage()->Add(m_DataNode);
    m_AblationSpheres.push_back(m_DataNode);

    // Add Center Points
    double finalRadius =
      m_AblationPlan->GetAblationZone(index)->radius *
      (1 - RadiusModellingUtils::calculateShrinkageOfARadiusDophi(m_AblationPlan->GetAblationZone(index)->radius));
    ;
    MITK_INFO << "Ablation Zone " << index << "[" << centerInWorldCoordinates[0] << ";" << centerInWorldCoordinates[1]
              << ";" << centerInWorldCoordinates[2] << "] / Radius: " << finalRadius;
    centerPoints->InsertPoint(centerInWorldCoordinates);
    zoneCenters.push_back(centerPoints);
  }
  m_AblationCenters = centerPoints;
  this->CreateContourBtwCenters(COG, centerPoints);
  m_AblationCentersNode = mitk::DataNode::New();
  m_AblationCentersNode->SetName("Ablation Centers");
  m_AblationCentersNode->SetData(centerPoints);
  // this->VisualizeMovedCenters(COG, centerPoints);
  this->GetDataStorage()->Add(m_AblationCentersNode);

  this->GetDataStorage()->Modified();
  this->RequestRenderWindowUpdate();
}

void QmitkAblationPlanningView::DeleteAllSpheres()
{
  for (int index = m_AblationZonesProcessed.size(); index > 0; --index)
  {
    QString name = QString("Kugel_%1").arg(index);
    mitk::DataNode::Pointer dataNode = this->GetDataStorage()->GetNamedNode(name.toStdString());
    if (dataNode.IsNotNull())
    {
      this->GetDataStorage()->Remove(dataNode);
    }
  }
  this->GetDataStorage()->Remove(m_AblationCentersNode);
  // this->GetDataStorage()->Remove(m_MovedCentersCOG);
  this->GetDataStorage()->Modified();
  for (mitk::DataNode::Pointer p : m_AblationSpheres)
  {
    this->GetDataStorage()->Remove(p);
  }
  m_AblationSpheres.clear();
  this->RequestRenderWindowUpdate();
}

void QmitkAblationPlanningView::DeleteContours()
{
  this->GetDataStorage()->Remove(m_NewNode);
  this->GetDataStorage()->Modified();
  for (mitk::DataNode::Pointer p : contourNodes)
  {
    this->GetDataStorage()->Remove(p);
    this->GetDataStorage()->Modified();
  }
  contourNodes.clear();
  this->RequestRenderWindowUpdate();
}

void QmitkAblationPlanningView::CalculateAblationStatistics()
{
  m_Controls.numberAblationVoluminaLabel->setText(QString("%1").arg(m_AblationPlan->GetNumberOfZones()));
  m_Controls.numberTumorVolumeLabel->setText(QString("%1").arg(m_AblationPlan->GetStatistics().tumorVolume));
  m_Controls.numberTumorAndMarginVolumeLabel->setText(
    QString("%1").arg(m_AblationPlan->GetStatistics().tumorAndSafetyMarginVolume));
  m_Controls.numberAblationVolumeLabel->setText(QString("%1").arg(m_AblationPlan->GetStatistics().totalAblationVolume));
  m_Controls.numberVolumeAblatedTwoAndMoreLabel->setText(
    QString("%1").arg(m_AblationPlan->GetStatistics().ablationVolumeAblatedMoreThanOneTime));
  m_Controls.numberOverlappingAblationZonesLabel->setText(
    QString("%1").arg(m_AblationPlan->GetStatistics().factorOverlappingAblationZones));
  m_Controls.numberFactorAblatedVolumeOutsideSafetyMarginLabel->setText(
    QString("%1").arg(m_AblationPlan->GetStatistics().factorAblatedVolumeOutsideSafetyMargin));
}

void QmitkAblationPlanningView::OnSegmentationComboBoxSelectionChanged(const mitk::DataNode *node)
{
  MITK_DEBUG << "OnSegmentationComboBoxSelectionChanged()";

  if (m_SegmentationImage.IsNotNull())
  {
    AblationUtils::ResetSafetyMargin(m_SegmentationImage, m_ImageDimension);
    AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension);
    this->DeleteAllSpheres();
    this->DeleteContours();

    m_TumorTissueSafetyMarginIndices.clear();
    m_AblationZones.clear();
    m_TempAblationZones.clear();
    m_AblationZonesProcessed.clear();
    m_TempAblationZonesProcessed.clear();

    this->CalculateAblationStatistics();
    mitk::RenderingManager::GetInstance()->Modified();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }

  if (node == nullptr)
  {
    this->UnsetSegmentationImageGeometry();
    m_SegmentationImage = nullptr;
    return;
  }

  mitk::DataNode *selectedSegmentation = m_Controls.segmentationComboBox->GetSelectedNode();
  if (selectedSegmentation == nullptr)
  {
    this->UnsetSegmentationImageGeometry();
    m_SegmentationImage = nullptr;
    return;
  }

  mitk::Image::Pointer segmentationImage = dynamic_cast<mitk::Image *>(selectedSegmentation->GetData());
  if (segmentationImage.IsNull())
  {
    MITK_WARN << "Failed to cast selected segmentation node to mitk::Image*";
    this->UnsetSegmentationImageGeometry();
    m_SegmentationImage = nullptr;
    return;
  }

  m_SegmentationImage = segmentationImage;
  this->SetSegmentationImageGeometryInformation(segmentationImage.GetPointer());
}

void QmitkAblationPlanningView::OnVisiblePropertyChanged()
{
  MITK_DEBUG << "OnVisiblePropertyChanged()";
}

void QmitkAblationPlanningView::OnBinaryPropertyChanged()
{
  MITK_DEBUG << "OnBinaryPropertyChanged()";
  mitk::DataStorage::SetOfObjects::ConstPointer segImages = m_Controls.segmentationComboBox->GetNodes();

  for (mitk::DataStorage::SetOfObjects::ConstIterator it = segImages->Begin(); it != segImages->End(); ++it)
  {
    const mitk::DataNode *node = it->Value();
    if (!m_IsASegmentationImagePredicate->CheckNode(node))
    {
      m_Controls.segmentationComboBox->RemoveNode(node);
      return;
    }
  }
}

void QmitkAblationPlanningView::OnCalculateSafetyMargin()
{
  if (m_SegmentationImage.IsNotNull() && m_Controls.safetyMarginSpinBox->value() > 0.0)
  {
    MITK_INFO << "Calculating safety margin is in progress...";
    // Reset the old calculated safety margin if the user calculated a safety margin before:
    AblationUtils::ResetSafetyMargin(m_SegmentationImage, m_ImageDimension);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          // Check, if von-Neumann neighbour elements are tumor tissue or not:
          if (!AblationUtils::CheckAllVonNeumannNeighbourPixelsAreTumorTissue(
                actualIndex, m_SegmentationImage, m_ImageDimension))
          {
            double margin = m_Controls.safetyMarginSpinBox->value();
            AblationUtils::CreateSafetyMarginInfluenceAreaOfPixel(
              actualIndex, m_SegmentationImage, margin, m_ImageDimension, m_ImageSpacing);
          }
        }
      }
    }
    MITK_INFO << "Finished calculating safety margin.";
    mitk::RenderingManager::GetInstance()->Modified();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
  else if (m_SegmentationImage.IsNotNull() && m_Controls.safetyMarginSpinBox->value() == 0.0)
  {
    AblationUtils::ResetSafetyMargin(m_SegmentationImage, m_ImageDimension);
    MITK_INFO << "Reset safety margin done.";
    mitk::RenderingManager::GetInstance()->Modified();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
  else
  {
    QMessageBox msgBox;
    msgBox.setText("Cannot calculate safety margin: No segmentation image was chosen.");
    msgBox.exec();
  }
}

void QmitkAblationPlanningView::CreateNodeForTumorCOG(mitk::Image::Pointer m_SegmentationImage)
{
  // get selected surface and name of surface
  selectedSurface = dynamic_cast<mitk::Surface *>(this->m_Controls.applySBox->GetSelectedNode()->GetData());
  std::string nameOfSelectedSurface = this->m_Controls.applySBox->GetSelectedNode()->GetName();

  mitk::PointSet::Pointer points = AblationUtils::CalculateCOGTargetPoints(selectedSurface);

  // create a new data node with targets
  m_NewNode = mitk::DataNode::New();
  m_NewNode->SetName(nameOfSelectedSurface + "_CenterOfGravity");
  m_NewNode->SetData(points);
  m_NewNode->SetProperty("color", mitk::ColorProperty::New(0, 255, 0));
  m_NewNode->SetProperty("pointsize", mitk::FloatProperty::New(3));
  m_NewNode->SetProperty("offset", mitk::FloatProperty::New(3));
  // add the new node to the data storage
  this->GetDataStorage()->Add(m_NewNode);
  this->GetDataStorage()->Modified();
}

void QmitkAblationPlanningView::CreateContourBtwCenters(mitk::PointSet::Pointer COG,
                                                        mitk::PointSet::Pointer zoneCenters)
{
  std::vector<mitk::PointSet::Pointer> pointSets(zoneCenters->GetSize());
  std::vector<mitk::DataNode::Pointer> pointSetNodes(zoneCenters->GetSize());
  // make new selection
  QList<mitk::DataNode::Pointer> selection;
  for (int index = 0; index < zoneCenters->GetSize(); ++index)
  {
    pointSets[index] = mitk::PointSet::New();
    pointSets[index]->InsertPoint(0, COG->GetPoint(0));
    pointSets[index]->InsertPoint(1, zoneCenters->GetPoint(index));
    pointSetNodes[index] = mitk::DataNode::New();
    pointSetNodes[index]->SetData(pointSets[index]);
    QString name = "Line COG to center" + QString::number(index + 1);
    pointSetNodes[index]->SetProperty("name", mitk::StringProperty::New(name.toStdString()));
    pointSetNodes[index]->SetProperty("show contour", mitk::BoolProperty::New(true));
    pointSetNodes[index]->SetProperty("show distances", mitk::BoolProperty::New(true));
    pointSetNodes[index]->SetProperty("contourcolor", mitk::ColorProperty::New(255, 0, 255));
    contourNodes.push_back(pointSetNodes[index]);
    this->GetDataStorage()->Add(contourNodes[index]);
    // on selection changed
    selection.push_back(pointSetNodes[index]);
    this->FireNodesSelected(selection);
    this->OnSelectionChanged(berry::IWorkbenchPart::Pointer(), selection);
  }
}

/* void QmitkAblationPlanningView::VisualizeMovedCenters(mitk::PointSet::Pointer COG,
                                                         mitk::PointSet::Pointer zoneCenters)
{
  m_MovedCentersCOG = mitk::DataNode::New();
  mitk::PointSet::Pointer newCoordinates = AblationUtils::GetCoordinatesBasedOnCOG(COG, zoneCenters);
  m_MovedCentersCOG->SetData(newCoordinates);
  m_MovedCentersCOG->SetName("COG-centers moved");
  m_MovedCentersCOG->SetProperty("pointsize", mitk::FloatProperty::New(5));
  this->GetDataStorage()->Add(m_MovedCentersCOG);
}*/

void QmitkAblationPlanningView::OnAblationStartingPointPushButtonClicked()
{
  if (m_SegmentationImage.IsNull())
  {
    MITK_WARN << "No segmentation image was chosen. Cannot set ablation starting point.";
    m_ManualAblationStartingPositionSet = false;
    return;
  }

  // Get the actual marked position of the crosshair:
  m_AblationStartingPositionInWorldCoordinates = this->GetRenderWindowPart()->GetSelectedPosition();

  // Calculate the index coordinates of the starting position:
  m_SegmentationImage->GetGeometry()->WorldToIndex(m_AblationStartingPositionInWorldCoordinates,
                                                   m_AblationStartingPositionIndexCoordinates);

  mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
  unsigned short pixelType = imagePixelWriter.GetPixelByIndex(m_AblationStartingPositionIndexCoordinates);

  MITK_DEBUG << "PixelType: " << pixelType;
  if (pixelType < 1)
  {
    m_ManualAblationStartingPositionSet = false;
    m_Controls.ablationStartingPointLabel->setText(
      "Position is not in the segmentation. Please choose a new starting position.");
    return;
  }

  m_ManualAblationStartingPositionSet = true;
  double x = m_AblationStartingPositionInWorldCoordinates[0];
  double y = m_AblationStartingPositionInWorldCoordinates[1];
  double z = m_AblationStartingPositionInWorldCoordinates[2];
  QString text = QString("Set Ablation Startingposition to: %1 | %2 | %3").arg(x).arg(y).arg(z);
  m_Controls.ablationStartingPointLabel->setText(text);
  MITK_DEBUG << "Set Ablation Startingposition to: " << m_AblationStartingPositionInWorldCoordinates;
  MITK_DEBUG << "Startingposition in Index: " << m_AblationStartingPositionIndexCoordinates;
  MITK_DEBUG << "Spacing: " << m_SegmentationImage->GetGeometry()->GetSpacing();
  // Get number of voxels in the three dimensions:
  MITK_DEBUG << "Dimension: " << m_ImageDimension[0] << " " << m_ImageDimension[1] << " " << m_ImageDimension[2];
}

void QmitkAblationPlanningView::OnCalculateAblationZonesPushButtonClicked()
{
  if (m_SegmentationImage.IsNull())
  {
    QMessageBox msgBox;
    msgBox.setText("Cannot calculate ablation zones. SegmentationImage is NULL.");
    msgBox.exec();
    return;
  }

  // check if a surface was created and the input for COG is valid
  if (this->m_Controls.applySBox->GetSelectedNode().IsNull())
  {
    MITK_WARN << "Error, no surface selected; Please create a smoothed surface";
    return;
  }

  m_AblationCalculationMade = true;

  // Reset earlier calculations:
  this->DeleteAllSpheres();
  this->DeleteContours();
  AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension);
  m_TumorTissueSafetyMarginIndices.clear();
  m_AblationZones.clear();
  m_TempAblationZones.clear();
  m_AblationZonesProcessed.clear();
  m_TempAblationZonesProcessed.clear();

  // Get some parameters from UI and set them for algorithm
  double minAblationRadius = RadiusModellingUtils::getPreAblationMinRadiusDophi();
  MITK_INFO << "MINNNNN " << minAblationRadius;
  double maxAblationRadius = RadiusModellingUtils::getPreAblationMaxRadiusDophi();
  MITK_INFO << "MAXXXXX " << maxAblationRadius;
  double toleranceNonAblatedVolume =
    (double)m_Controls.toleranceNonAblatedTumorSafetyMarginVolumeSpinBox->value() / 100.0;
  m_PlanningAlgo->SetAdjustableParameters(m_Controls.repititionsCalculatingAblationZonesSpinBox->value(),
                                          maxAblationRadius,
                                          minAblationRadius,
                                          toleranceNonAblatedVolume);
  m_PlanLogger->SetFileName(m_Controls.m_LoggingFileName->text().toStdString());

  // Set segmentation image for algorithm
  m_PlanningAlgo->SetSegmentationData(m_SegmentationImage, m_ImageDimension, m_ImageSpacing);

  // Set safety margin for algorithm
  m_PlanningAlgo->SetSafetyMargin(m_TumorTissueSafetyMarginIndices);

  // Compute planning
  m_PlanningAlgo->ComputePlanning();
  MITK_INFO << "Finished calculating ablation zones!";

  // Get final proposal and visualize it!
  mitk::AblationPlan::Pointer finalProposal = m_PlanningAlgo->GetAblationPlan();
  m_Controls.numberAblationVoluminaLabel->setText(QString::number(finalProposal->GetNumberOfZones()));
  mitk::RenderingManager::GetInstance()->Modified();
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  double notAblated = AblationUtils::CheckImageForNonAblatedTissueInPercentage(finalProposal->GetSegmentationImage(),
                                                                               finalProposal->GetImageDimension());
  if (notAblated > 0)
  {
    MITK_WARN << "There is still non ablated tumor tissue (" << notAblated << " percent).";
  }
  m_AblationPlan = finalProposal;
  this->CreateNodeForTumorCOG(m_SegmentationImage);
  this->CreateSpheresOfAblationVolumes();
  this->CalculateAblationStatistics();

  // To get the tumor COG for the log file
  mitk::PointSet::Pointer COG = AblationUtils::CalculateCOGTargetPoints(selectedSurface);

  // To get the coordinates of the ablation centers relativ to the tumor COG
  mitk::PointSet::Pointer newCoordinates = AblationUtils::GetCoordinatesBasedOnCOG(COG, m_AblationCenters);

  // logging of results
  if (this->m_Controls.LoggingActivated->isChecked())
  {
    std::stringstream caseName;
    caseName << this->m_Controls.segmentationComboBox->GetSelectedNode()->GetName();
    std::time_t t = std::time(0);
    std::tm *now = std::localtime(&t);
    caseName << "y" << (now->tm_year + 1900) << "m" << now->tm_mon + 1 << "d" << now->tm_mday << "h" << now->tm_hour
             << "m" << now->tm_min << "s" << now->tm_sec;
    mitk::AblationPlanningLogging::AblationPlanningParameterSet params;
    params.maxRadius = RadiusModellingUtils::getMaxRadiusDophi();
    params.minRadius = RadiusModellingUtils::getMinRadiusDophi();
    params.iterations = m_Controls.repititionsCalculatingAblationZonesSpinBox->value();
    params.safetyMargin = m_Controls.safetyMarginSpinBox->value();
    // params.tissueShrinking = 0.2;
    params.toleranceNonAblatedVolume = m_Controls.toleranceNonAblatedTumorSafetyMarginVolumeSpinBox->value();
    params.COGravity = COG;
    params.relativeCoordinates = newCoordinates;
    m_PlanLogger->WriteHeader();
    m_PlanLogger->WriteDataSet(
      finalProposal, m_Controls.segmentationComboBox->GetSelectedNode(), params, caseName.str());
    m_PlanLogger->WriteScene(this->GetDataStorage(), caseName.str());
    MITK_INFO << "Logged all results to file " << m_Controls.m_LoggingFileName->text().toStdString() << " under name "
              << caseName.str();
  }
}

void QmitkAblationPlanningView::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  m_Controls.segmentationComboBox->SetDataStorage(GetDataStorage());
  m_Controls.segmentationComboBox->SetPredicate(m_IsASegmentationImagePredicate);
  // connecting the surface from data manager with surface combo box
  m_Controls.applySBox->SetDataStorage(this->GetDataStorage());
  m_Controls.applySBox->SetAutoSelectNewItems(true);
  m_Controls.applySBox->SetPredicate(mitk::NodePredicateDataType::New("Surface"));
  // create signal/slot connections
  connect(m_Controls.segmentationComboBox,
          SIGNAL(OnSelectionChanged(const mitk::DataNode *)),
          this,
          SLOT(OnSegmentationComboBoxSelectionChanged(const mitk::DataNode *)));
  connect(m_Controls.ablationStartingPointPushButton,
          SIGNAL(clicked()),
          this,
          SLOT(OnAblationStartingPointPushButtonClicked()));
  connect(m_Controls.calculateAblationZonesPushButton,
          SIGNAL(clicked()),
          this,
          SLOT(OnCalculateAblationZonesPushButtonClicked()));
  connect(m_Controls.calculateSafetyMarginPushButton, SIGNAL(clicked()), this, SLOT(OnCalculateSafetyMargin()));
  connect(m_Controls.m_ChooseFile, SIGNAL(clicked()), this, SLOT(OnChooseFileClicked()));

  m_Controls.m_LoggingFileName->setText(QDir::toNativeSeparators(((QDir)QDir::homePath()).absolutePath()) +
                                        QDir::separator() + "output.csv");
  mitk::DataStorage::SetOfObjects::ConstPointer segmentationImages =
    GetDataStorage()->GetSubset(m_IsASegmentationImagePredicate);
  if (!segmentationImages->empty())
  {
    OnSelectionChanged(*segmentationImages->begin());
  }

  // set callback function for already existing nodes (segmentations)
  mitk::DataStorage::SetOfObjects::ConstPointer allImages = GetDataStorage()->GetSubset(m_IsOfTypeImagePredicate);
  for (mitk::DataStorage::SetOfObjects::const_iterator iter = allImages->begin(); iter != allImages->end(); ++iter)
  {
    mitk::DataNode *node = *iter;
    itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command =
      itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
    command->SetCallbackFunction(this, &QmitkAblationPlanningView::OnVisiblePropertyChanged);
    m_WorkingDataObserverTags.insert(std::pair<mitk::DataNode *, unsigned long>(
      node, node->GetProperty("visible")->AddObserver(itk::ModifiedEvent(), command)));

    itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command2 =
      itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
    command2->SetCallbackFunction(this, &QmitkAblationPlanningView::OnBinaryPropertyChanged);
    m_BinaryPropertyObserverTags.insert(std::pair<mitk::DataNode *, unsigned long>(
      node, node->GetProperty("binary")->AddObserver(itk::ModifiedEvent(), command2)));
  }
}

void QmitkAblationPlanningView::OnChooseFileClicked()
{
  QDir currentPath = QFileInfo(m_Controls.m_LoggingFileName->text()).dir();

  // if no path was selected (QDir would select current working dir then) or the
  // selected path does not exist -> use home directory
  if (currentPath == QDir() || !currentPath.exists())
  {
    currentPath = QDir(QDir::homePath());
  }

  QString filename =
    QFileDialog::getSaveFileName(nullptr, tr("Choose Logging File"), currentPath.absolutePath(), "*.csv");
  if (filename == "")
    return;
  m_Controls.m_LoggingFileName->setText(filename);
}
