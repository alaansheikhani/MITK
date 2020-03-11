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
#include <QFileDialog>
#include <QTextStream>

// mitk
#include "mitkAblationUtils.h"

#include <mitkImage.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkLabelSetImage.h>
#include <mitkPointSet.h>
#include "mitkProperties.h"
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkSceneIO.h>

#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <mitkSurface.h>
#include <chrono>

const std::string QmitkAblationPlanningView::VIEW_ID = "org.mitk.views.ablationplanning";
const static short ABLATION_VALUE = 2;
const static short TUMOR_NOT_YET_ABLATED = 1;
const static short NO_TUMOR_ISSUE = 0;
const static short SAFETY_MARGIN = 256;
const static unsigned short BIT_OPERATION_ELIMINATE_TUMOR_SAFETY_MARGIN = 65278; // = 11111110 11111110

//=====================Konstruktor/Destruktor===================================
QmitkAblationPlanningView::QmitkAblationPlanningView()
  :
  m_MouseCursorSet(false),
  m_DataSelectionChanged(false),
  m_AblationStartingPositionInWorldCoordinates(),
  m_AblationStartingPositionIndexCoordinates(),
  m_ManualAblationStartingPositionSet(false),
  m_AblationRadius(15.0),
  m_AblationCalculationMade(false),
  m_AblationCentersNode(mitk::DataNode::New())
{
  this->UnsetSegmentationImageGeometry();

  mitk::TNodePredicateDataType<mitk::Image>::Pointer isImage = mitk::TNodePredicateDataType<mitk::Image>::New();
  //mitk::NodePredicateDataType::Pointer isDwi = mitk::NodePredicateDataType::New("DiffusionImage");
  //mitk::NodePredicateDataType::Pointer isDti = mitk::NodePredicateDataType::New("TensorImage");
  //mitk::NodePredicateDataType::Pointer isOdf = mitk::NodePredicateDataType::New("OdfImage");
  auto isSegmentation = mitk::NodePredicateDataType::New("Segment");

  mitk::NodePredicateOr::Pointer validImages = mitk::NodePredicateOr::New();
  validImages->AddPredicate(mitk::NodePredicateAnd::New(isImage, mitk::NodePredicateNot::New(isSegmentation)));
  //validImages->AddPredicate(isDwi);
  //validImages->AddPredicate(isDti);
  //validImages->AddPredicate(isOdf);

  mitk::NodePredicateNot::Pointer isNotAHelperObject = mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true)));

  m_IsOfTypeImagePredicate = mitk::NodePredicateAnd::New(validImages, isNotAHelperObject);

  mitk::NodePredicateProperty::Pointer isBinaryPredicate = mitk::NodePredicateProperty::New("binary", mitk::BoolProperty::New(true));
  mitk::NodePredicateNot::Pointer isNotBinaryPredicate = mitk::NodePredicateNot::New(isBinaryPredicate);

  mitk::NodePredicateAnd::Pointer isABinaryImagePredicate = mitk::NodePredicateAnd::New(m_IsOfTypeImagePredicate, isBinaryPredicate);
  mitk::NodePredicateAnd::Pointer isNotABinaryImagePredicate = mitk::NodePredicateAnd::New(m_IsOfTypeImagePredicate, isNotBinaryPredicate);

  m_IsASegmentationImagePredicate = mitk::NodePredicateOr::New(isABinaryImagePredicate, mitk::TNodePredicateDataType<mitk::LabelSetImage>::New());
  m_IsAPatientImagePredicate = mitk::NodePredicateAnd::New(isNotABinaryImagePredicate, mitk::NodePredicateNot::New(mitk::TNodePredicateDataType<mitk::LabelSetImage>::New()));

  // to initialize with the first image, if nothing was selected yet...
  if (this->GetDataStorage()->GetNode(m_IsASegmentationImagePredicate) != nullptr) {
    m_SegmentationImage = dynamic_cast<mitk::Image*>(this->GetDataStorage()->GetNode(m_IsASegmentationImagePredicate)->GetData());
    SetSegmentationImageGeometryInformation(m_SegmentationImage);
    OnSelectionChanged(this->GetDataStorage()->GetNode(m_IsASegmentationImagePredicate));
  }
}

QmitkAblationPlanningView::~QmitkAblationPlanningView()
{
  // removing all observers
  for (NodeTagMapType::iterator dataIter = m_WorkingDataObserverTags.begin(); dataIter != m_WorkingDataObserverTags.end(); ++dataIter)
  {
    (*dataIter).first->GetProperty("visible")->RemoveObserver((*dataIter).second);
  }
  m_WorkingDataObserverTags.clear();

  for (NodeTagMapType::iterator dataIter = m_BinaryPropertyObserverTags.begin(); dataIter != m_BinaryPropertyObserverTags.end(); ++dataIter)
  {
    (*dataIter).first->GetProperty("binary")->RemoveObserver((*dataIter).second);
  }
  m_BinaryPropertyObserverTags.clear();

  //mitk::RenderingManager::GetInstance()->RemoveObserver(m_RenderingManagerObserverTag);

}

void QmitkAblationPlanningView::OnSelectionChanged(mitk::DataNode * node)
{
  MITK_DEBUG << "OnSelectionChanged()";
  berry::IWorkbenchPart::Pointer nullPart;
  QList<mitk::DataNode::Pointer> nodes;
  nodes.push_back(node);
  this->OnSelectionChanged(nullPart, nodes);
}

void QmitkAblationPlanningView::OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes)
{
  MITK_DEBUG << "OnSelectionChanged()";

  if (nodes.size() == 1)
  {
    mitk::DataNode::Pointer selectedNode = nodes.at(0);
    if (selectedNode.IsNull())
    {
      return;
    }

    mitk::Image::Pointer selectedImage = dynamic_cast<mitk::Image*>(selectedNode->GetData());
    if (selectedImage.IsNull())
    {
      return;
    }

    if (m_IsASegmentationImagePredicate->CheckNode(selectedNode))
    {
      // if a segmentation is selected find a possible patient image
      mitk::DataStorage::SetOfObjects::ConstPointer sources = GetDataStorage()->GetSources(selectedNode, m_IsAPatientImagePredicate);
      mitk::DataNode::Pointer refNode;
      if (sources->Size() != 0)
      {
        // found one or more sources - use the first one
        refNode = sources->ElementAt(0);
        refNode->SetVisibility(true);
        selectedNode->SetVisibility(true);
      }
      //mitk::RenderingManager::GetInstance()->InitializeViews(selectedNode->GetData()->GetTimeGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true);
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

void QmitkAblationPlanningView::SetSegmentationImageGeometryInformation(mitk::Image* image)
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

void QmitkAblationPlanningView::NodeRemoved(const mitk::DataNode * node)
{
  MITK_DEBUG << "NodeRemoved()";
  if (m_IsASegmentationImagePredicate->CheckNode(node))
  {
    //First of all remove all possible contour markers of the segmentation
    mitk::DataStorage::SetOfObjects::ConstPointer allContourMarkers = this->GetDataStorage()->GetDerivations(node, mitk::NodePredicateProperty::New("isContourMarker", mitk::BoolProperty::New(true)));

    for (mitk::DataStorage::SetOfObjects::ConstIterator it = allContourMarkers->Begin(); it != allContourMarkers->End(); ++it)
    {
      std::string nodeName = node->GetName();
      unsigned int t = nodeName.find_last_of(" ");
      unsigned int id = atof(nodeName.substr(t + 1).c_str()) - 1;

      this->GetDataStorage()->Remove(it->Value());
    }

    mitk::Image* image = dynamic_cast<mitk::Image*>(node->GetData());
  }
  mitk::DataNode* tempNode = const_cast<mitk::DataNode*>(node);
  //Since the binary property could be changed during runtime by the user
  if (m_IsOfTypeImagePredicate->CheckNode(node))
  {
    node->GetProperty("visible")->RemoveObserver(m_WorkingDataObserverTags[tempNode]);
    m_WorkingDataObserverTags.erase(tempNode);
    node->GetProperty("binary")->RemoveObserver(m_BinaryPropertyObserverTags[tempNode]);
    m_BinaryPropertyObserverTags.erase(tempNode);
  }
}

void QmitkAblationPlanningView::NodeAdded(const mitk::DataNode * node)
{
  MITK_DEBUG << "NodeAdded()";
  if (!m_IsOfTypeImagePredicate->CheckNode(node))
  {
    return;
  }

  itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command = itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
  command->SetCallbackFunction(this, &QmitkAblationPlanningView::OnVisiblePropertyChanged);
  m_WorkingDataObserverTags.insert(std::pair<mitk::DataNode*, unsigned long>(const_cast<mitk::DataNode*>(node), node->GetProperty("visible")->AddObserver(itk::ModifiedEvent(), command)));

  itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command2 = itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
  command2->SetCallbackFunction(this, &QmitkAblationPlanningView::OnBinaryPropertyChanged);
  m_BinaryPropertyObserverTags.insert(std::pair<mitk::DataNode*, unsigned long>(const_cast<mitk::DataNode*>(node), node->GetProperty("binary")->AddObserver(itk::ModifiedEvent(), command2)));
}

bool QmitkAblationPlanningView::CheckForSameGeometry(const mitk::DataNode *node1, const mitk::DataNode *node2) const
{
  bool isSameGeometry(true);

  mitk::Image* image1 = dynamic_cast<mitk::Image*>(node1->GetData());
  mitk::Image* image2 = dynamic_cast<mitk::Image*>(node2->GetData());
  if (image1 && image2)
  {
    mitk::BaseGeometry* geo1 = image1->GetGeometry();
    mitk::BaseGeometry* geo2 = image2->GetGeometry();

    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetOrigin(), geo2->GetOrigin());
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(0), geo2->GetExtent(0));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(1), geo2->GetExtent(1));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetExtent(2), geo2->GetExtent(2));
    isSameGeometry = isSameGeometry && mitk::Equal(geo1->GetSpacing(), geo2->GetSpacing());
    isSameGeometry = isSameGeometry && mitk::MatrixEqualElementWise(geo1->GetIndexToWorldTransform()->GetMatrix(), geo2->GetIndexToWorldTransform()->GetMatrix());

    return isSameGeometry;
  }
  else
  {
    return false;
  }
}

void QmitkAblationPlanningView::CopyTemporaryAblationZoneDistribution()
{
  if( m_AblationZoneCentersProcessed.size() == 0 ||
      m_AblationZoneCentersProcessed.size() > m_TempAblationZoneCentersProcessed.size())
  {
    MITK_INFO << "Reduced the number of ablation zones from: "
              << m_AblationZoneCentersProcessed.size()
              << " to: "
              << m_TempAblationZoneCentersProcessed.size()
              << "========================================";

    m_AblationZoneCenters.clear();
    m_AblationZoneCentersProcessed.clear();
    m_AblationZoneCenters = m_TempAblationZoneCenters;
    m_AblationZoneCentersProcessed = m_TempAblationZoneCentersProcessed;
    m_TempAblationZoneCenters.clear();
    m_TempAblationZoneCentersProcessed.clear();
    m_AblationStartingPositionIndexCoordinates = m_TempAblationStartingPositionIndexCoordinates;
    m_AblationStartingPositionInWorldCoordinates = m_TempAblationStartingPositionInWorldCoordinates;
  }
  else
  {
    m_TempAblationZoneCenters.clear();
    m_TempAblationZoneCentersProcessed.clear();
  }
}

void QmitkAblationPlanningView::CreateSpheresOfAblationVolumes()
{
  mitk::PointSet::Pointer centerPoints = mitk::PointSet::New();
  for( int index = 0; index < m_AblationZoneCentersProcessed.size(); ++index)
  {
    mitk::DataNode::Pointer m_DataNode = mitk::DataNode::New();

    mitk::Surface::Pointer mySphere = mitk::Surface::New();

    vtkSmartPointer<vtkSphereSource> vtkSphere = vtkSmartPointer<vtkSphereSource>::New();
    vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    mitk::Point3D centerInWorldCoordinates;

    m_SegmentationImage->GetGeometry()->IndexToWorld(
                                          m_AblationZoneCentersProcessed.at(index),
                                          centerInWorldCoordinates );

    //Center
    vtkSphere->SetRadius(m_AblationRadius);
    vtkSphere->SetPhiResolution(120);
    vtkSphere->SetThetaResolution(120);
    vtkSphere->SetCenter(centerInWorldCoordinates[0], centerInWorldCoordinates[1], centerInWorldCoordinates[2]);
    vtkSphere->Update();

    appendPolyData->AddInputData(vtkSphere->GetOutput());

    mySphere->SetVtkPolyData(vtkSphere->GetOutput());

    //Add Node
    m_DataNode->SetData(mySphere);
    QString name = QString("Kugel_%1").arg(index + 1);
    m_DataNode->SetName(name.toStdString());
    this->GetDataStorage()->Add(m_DataNode);

    //Add Center Points
    MITK_INFO << "Ablation Zone " << index << "[" << centerInWorldCoordinates[0] << ";" << centerInWorldCoordinates[1]
              << ";" << centerInWorldCoordinates[2] << "]";
    centerPoints->InsertPoint(centerInWorldCoordinates);
  }

  m_AblationCentersNode = mitk::DataNode::New();
  m_AblationCentersNode->SetName("Ablation Centers");
  m_AblationCentersNode->SetData(centerPoints);
  this->GetDataStorage()->Add(m_AblationCentersNode);

  this->GetDataStorage()->Modified();
  this->RequestRenderWindowUpdate();
}

void QmitkAblationPlanningView::DeleteAllSpheres()
{
  for( int index = m_AblationZoneCentersProcessed.size(); index > 0; --index )
  {
    QString name = QString("Kugel_%1").arg(index);

    mitk::DataNode::Pointer dataNode = this->GetDataStorage()->GetNamedNode(name.toStdString());
    if( dataNode.IsNotNull() )
    {
      this->GetDataStorage()->Remove(dataNode);
    }
  }
  this->GetDataStorage()->Remove(m_AblationCentersNode);
  this->GetDataStorage()->Modified();
  this->RequestRenderWindowUpdate();
}

void QmitkAblationPlanningView::FillComboBoxAblationZones()
{
  m_Controls.ablationZonesComboBox->clear();

  for (int index = 0; index < m_AblationZoneCentersProcessed.size(); ++index)
  {
    QString name = QString("Kugel_%1").arg(index + 1);
    m_Controls.ablationZonesComboBox->addItem(name);
  }
}

void QmitkAblationPlanningView::CalculateAblationStatistics()
{
  m_Stats.tumorVolume =
    AblationUtils::CalculateTumorVolume(m_SegmentationImage, m_ImageSpacing, m_TumorTissueSafetyMarginIndices);
  m_Stats.safetyMarginVolume =
    AblationUtils::CalculateSafetyMarginVolume(m_SegmentationImage, m_ImageSpacing, m_TumorTissueSafetyMarginIndices);
  m_Stats.tumorAndSafetyMarginVolume = m_Stats.tumorVolume + m_Stats.safetyMarginVolume;
  m_Stats.totalAblationVolume =
    AblationUtils::CalculateTotalAblationVolume(m_SegmentationImage, m_ImageSpacing, m_ImageDimension);
  m_Stats.ablationVolumeAblatedMoreThanOneTime =
    AblationUtils::CalculateAblationVolumeAblatedMoreThanOneTime(m_SegmentationImage, m_ImageSpacing, m_ImageDimension);

  m_Stats.factorOverlappingAblationZones =
    ((double)m_Stats.ablationVolumeAblatedMoreThanOneTime / m_Stats.totalAblationVolume) * 100;

  m_Stats.factorAblatedVolumeOutsideSafetyMargin =
    ((m_Stats.totalAblationVolume - m_Stats.tumorAndSafetyMarginVolume) / (double)m_Stats.totalAblationVolume) * 100;

  m_Stats.numberOfZones = m_AblationZoneCentersProcessed.size();
  m_Controls.numberAblationVoluminaLabel
    ->setText(QString("%1").arg(m_AblationZoneCentersProcessed.size()));
  m_Controls.numberTumorVolumeLabel->setText(QString("%1").arg(m_Stats.tumorVolume));
  m_Controls.numberTumorAndMarginVolumeLabel->setText(QString("%1").arg(m_Stats.tumorAndSafetyMarginVolume));
  m_Controls.numberAblationVolumeLabel->setText(QString("%1").arg(m_Stats.totalAblationVolume));
  m_Controls.numberVolumeAblatedTwoAndMoreLabel->setText(
    QString("%1").arg(m_Stats.ablationVolumeAblatedMoreThanOneTime));
  m_Controls.numberOverlappingAblationZonesLabel->setText(QString("%1").arg(m_Stats.factorOverlappingAblationZones));
  m_Controls.numberFactorAblatedVolumeOutsideSafetyMarginLabel->setText(
    QString("%1").arg(m_Stats.factorAblatedVolumeOutsideSafetyMargin));
}

void QmitkAblationPlanningView::OnSegmentationComboBoxSelectionChanged(const mitk::DataNode* node)
{
  MITK_DEBUG << "OnSegmentationComboBoxSelectionChanged()";

  if( m_SegmentationImage.IsNotNull() )
  {
    AblationUtils::ResetSafetyMargin(m_SegmentationImage, m_ImageDimension);
    AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension);
    this->DeleteAllSpheres();

    m_TumorTissueSafetyMarginIndices.clear();
    m_AblationZoneCenters.clear();
    m_TempAblationZoneCenters.clear();
    m_AblationZoneCentersProcessed.clear();
    m_TempAblationZoneCentersProcessed.clear();

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

  mitk::DataNode* selectedSegmentation = m_Controls.segmentationComboBox->GetSelectedNode();
  if (selectedSegmentation == nullptr)
  {
    this->UnsetSegmentationImageGeometry();
    m_SegmentationImage = nullptr;
    return;
  }

  mitk::Image::Pointer segmentationImage = dynamic_cast<mitk::Image*>(selectedSegmentation->GetData());
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
    const mitk::DataNode* node = it->Value();
    if (!m_IsASegmentationImagePredicate->CheckNode(node))
    {
      m_Controls.segmentationComboBox->RemoveNode(node);
      return;
    }
  }
}

void QmitkAblationPlanningView::OnMergeTumorSafetyMargin() {
  MITK_INFO << "Merging tumor segmentation an safety margin";
  AblationUtils::MergeSegmentationAndSecurityMargin(m_SegmentationImage, m_ImageDimension);
}

void QmitkAblationPlanningView::OnCalculateSafetyMargin()
{
  if(m_SegmentationImage.IsNotNull() && m_Controls.safetyMarginSpinBox->value() > 0.0 )
  {
    MITK_INFO << "Calculating safety margin is in progress...";
    //Reset the old calculated safety margin if the user calculated a safety margin before:
    AblationUtils::ResetSafetyMargin(m_SegmentationImage, m_ImageDimension);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          // Check, if von-Neumann neighbour elements are tumor tissue or not:
          if( !AblationUtils::CheckAllVonNeumannNeighbourPixelsAreTumorTissue(actualIndex, m_SegmentationImage, m_ImageDimension) )
          {
            double margin = m_Controls.safetyMarginSpinBox->value();
            AblationUtils::CreateSafetyMarginInfluenceAreaOfPixel(actualIndex, m_SegmentationImage, margin, m_ImageDimension, m_ImageSpacing);
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

void QmitkAblationPlanningView::OnAblationStartingPointPushButtonClicked()
{
/*
  mitk::DataNode* selectedSegmentation = m_Controls.segmentationComboBox->GetSelectedNode();
  if( selectedSegmentation == nullptr )
  {
    m_ManualAblationStartingPositionSet = false;
    return;
  }

  m_SegmentationImage = dynamic_cast<mitk::Image*>(selectedSegmentation->GetData());
  */

  if (m_SegmentationImage.IsNull())
  {
    MITK_WARN << "No segmentation image was chosen. Cannot set ablation starting point.";
    m_ManualAblationStartingPositionSet = false;
    return;
  }

  //Get the actual marked position of the crosshair:
  m_AblationStartingPositionInWorldCoordinates = this->GetRenderWindowPart()->GetSelectedPosition();

  //Calculate the index coordinates of the starting position:
  m_SegmentationImage->GetGeometry()->WorldToIndex( m_AblationStartingPositionInWorldCoordinates,
                                                    m_AblationStartingPositionIndexCoordinates);

  mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
  unsigned short pixelType = imagePixelWriter.GetPixelByIndex(m_AblationStartingPositionIndexCoordinates);

  MITK_DEBUG << "PixelType: " << pixelType;
  if (pixelType < 1)
  {
    m_ManualAblationStartingPositionSet = false;
    m_Controls.ablationStartingPointLabel->setText("Position is not in the segmentation. Please choose a new starting position.");
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
  //Get number of voxels in the three dimensions:
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

  m_AblationCalculationMade = true;
  m_Controls.refreshCalculationsPushButton->setVisible(false);

  //Reset earlier calculations:
  this->DeleteAllSpheres();
  AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension);
  m_TumorTissueSafetyMarginIndices.clear();
  m_AblationZoneCenters.clear();
  m_TempAblationZoneCenters.clear();
  m_AblationZoneCentersProcessed.clear();
  m_TempAblationZoneCentersProcessed.clear();

  std::chrono::milliseconds start = std::chrono::duration_cast< std::chrono::milliseconds >(
    std::chrono::system_clock::now().time_since_epoch());

  //==================== Find ablation proposal by iteratively create and test random propsals ==========================================
  AblationUtils::FillVectorContainingIndicesOfTumorTissueSafetyMargin(m_SegmentationImage, m_ImageDimension, m_TumorTissueSafetyMarginIndices);

  //Start of for-loop (main iterative loop, each iteration = one proposal):
  for( int iteration = 1; iteration <= m_Controls.repititionsCalculatingAblationZonesSpinBox->value(); ++iteration )
  {
    MITK_INFO << "Iteration: " << iteration;
    if (!m_ManualAblationStartingPositionSet || iteration > 1 )
    {
      QString position = AblationUtils::FindAblationStartingPosition(m_SegmentationImage, m_TumorTissueSafetyMarginIndices, m_AblationRadius, m_TempAblationStartingPositionIndexCoordinates, m_TempAblationStartingPositionInWorldCoordinates, m_ImageDimension, m_ImageSpacing);
      m_Controls.ablationStartingPointLabel->setText(position);
    }

    //------------ Grid calculation model: ---------------
    if( m_Controls.gridModelRadioButton->isChecked() )
    {
      AblationUtils::CalculateAblationVolume(m_TempAblationStartingPositionIndexCoordinates, m_SegmentationImage, m_AblationRadius, m_ImageSpacing, m_ImageDimension, m_TempAblationZoneCenters);
      AblationUtils::ProcessDirectNeighbourAblationZones(m_TempAblationStartingPositionIndexCoordinates, m_SegmentationImage, m_ImageSpacing, m_ImageDimension, m_AblationRadius, m_TempAblationZoneCentersProcessed, m_TempAblationZoneCenters);
      while( m_TempAblationZoneCenters.size() != m_TempAblationZoneCentersProcessed.size() )
      {
        MITK_DEBUG << "Size1: " << m_TempAblationZoneCenters.size()
                   << " Size2: " << m_TempAblationZoneCentersProcessed.size();
        for( int index = 0; index < m_TempAblationZoneCenters.size(); ++index )
        {
          if (!AblationUtils::IsAblationZoneAlreadyProcessed(m_TempAblationZoneCenters.at(index), m_TempAblationZoneCentersProcessed))
          {
            AblationUtils::ProcessDirectNeighbourAblationZones(m_TempAblationZoneCenters.at(index), m_SegmentationImage, m_ImageSpacing, m_ImageDimension, m_AblationRadius, m_TempAblationZoneCentersProcessed, m_TempAblationZoneCenters);
            break;
          }
        }
      }
    }
    //------------ Random distribution calculation model: ---------------
    else if( m_Controls.randomDistributionRadioButton->isChecked() )
    {
      double size = m_TumorTissueSafetyMarginIndices.size();
      std::vector<itk::Index<3>> indices = m_TumorTissueSafetyMarginIndices;
      AblationUtils::CalculateAblationVolume(m_TempAblationStartingPositionIndexCoordinates, m_SegmentationImage, m_AblationRadius, m_ImageSpacing, m_ImageDimension, m_TempAblationZoneCenters);
      AblationUtils::RemoveAblatedPixelsFromGivenVector(
        m_TempAblationStartingPositionIndexCoordinates, indices, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
      m_TempAblationZoneCentersProcessed.push_back(m_TempAblationStartingPositionIndexCoordinates);

      while( indices.size() > 0 &&
             (double)(indices.size() / size) >
             ((double)m_Controls.toleranceNonAblatedTumorSafetyMarginVolumeSpinBox->value() / 100) )
      {
        itk::Index<3> newAblationCenter = AblationUtils::SearchNextAblationCenter(indices, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
        AblationUtils::CalculateAblationVolume(newAblationCenter, m_SegmentationImage, m_AblationRadius, m_ImageSpacing, m_ImageDimension, m_TempAblationZoneCenters);
        AblationUtils::RemoveAblatedPixelsFromGivenVector(newAblationCenter, indices, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
        m_TempAblationZoneCentersProcessed.push_back(newAblationCenter);
      }
    }
    //------------ End calculation models -------------------------------

    AblationUtils::DetectNotNeededAblationVolume(m_TempAblationZoneCentersProcessed, m_TempAblationZoneCenters, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
    MITK_INFO << "Total number of ablation zones: " << m_TempAblationZoneCentersProcessed.size();
    this->CopyTemporaryAblationZoneDistribution();
    AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension);
  } //End of for loop

  //==================== Optimization of final proposal ==================================================

  //Check if ablation zones are too far outside the tumor, if yes move them towards the center
  for (int index = 0; index < m_AblationZoneCentersProcessed.size(); ++index)
  {
    double ratio =
      AblationUtils::CalculateRatioAblatedTissueOutsideTumorToAblatedTissueInsideTumor(
        m_AblationZoneCentersProcessed.at(index), m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
    MITK_WARN << "RATIO: " << ratio;
    if (ratio > 0.3)
    {
      AblationUtils::MoveCenterOfAblationZone(m_AblationZoneCentersProcessed.at(index), m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
    }
  }

  //Check, if ablation zones have a too short distance between each other, if yes they can be removed
  for(int index = 0; index < m_AblationZoneCentersProcessed.size(); ++index)
  {
    std::vector<int> indexToRemove;
    itk::Index<3> actualIndex = m_AblationZoneCentersProcessed.at(index);
    for (int counter = 0; counter < m_AblationZoneCentersProcessed.size(); ++counter)
    {
      if (counter == index)
      {
        continue;
      }
      itk::Index<3> indexToProof = m_AblationZoneCentersProcessed.at(counter);
      double distance = AblationUtils::CalculateScalarDistance(actualIndex, indexToProof, m_ImageSpacing);
      if( distance <= 0.5 * m_AblationRadius )
      {
        indexToRemove.push_back(counter);
      }
    }
    for (int position = indexToRemove.size() - 1; position >= 0; --position)
    {
      std::vector<itk::Index<3>>::iterator it = m_AblationZoneCentersProcessed.begin();
      m_AblationZoneCentersProcessed.erase(it + indexToRemove.at(position));
      std::vector<itk::Index<3>>::iterator it2 = m_AblationZoneCenters.begin();
      m_AblationZoneCenters.erase(it2 + indexToRemove.at(position));
      MITK_DEBUG << "Removed Ablation zone at index position: " << indexToRemove.at(position);
      index = -1;
    }

  }



  for (int index = 0; index < m_AblationZoneCentersProcessed.size(); ++index)
  {
    AblationUtils::CalculateAblationVolume(m_AblationZoneCentersProcessed.at(index), m_SegmentationImage, m_AblationRadius, m_ImageSpacing, m_ImageDimension, m_TempAblationZoneCenters);
  }

  AblationUtils::DetectNotNeededAblationVolume(m_AblationZoneCentersProcessed, m_AblationZoneCenters, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);


  //============
  if(m_Controls.toleranceNonAblatedTumorSafetyMarginVolumeSpinBox->value() == 0)
  {
    std::vector<itk::Index<3>> onlyTumorIndices = AblationUtils::FillVectorContainingIndicesOfTumorTissueOnly(m_SegmentationImage, m_ImageDimension);

    while (onlyTumorIndices.size() > 0)
    {
      itk::Index<3> newAblationCenter = AblationUtils::SearchNextAblationCenter(onlyTumorIndices, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
      AblationUtils::MoveCenterOfAblationZone(newAblationCenter, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
      AblationUtils::CalculateAblationVolume(newAblationCenter, m_SegmentationImage, m_AblationRadius, m_ImageSpacing, m_ImageDimension);
      AblationUtils::RemoveAblatedPixelsFromGivenVector(newAblationCenter, onlyTumorIndices, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
      m_AblationZoneCentersProcessed.push_back(newAblationCenter);
      m_AblationZoneCenters.push_back(newAblationCenter);
    }
    AblationUtils::DetectNotNeededAblationVolume(m_AblationZoneCentersProcessed, m_AblationZoneCenters, m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);
  }
  //============

  MITK_INFO << "Finished calculating ablation zones!";

  m_Controls.numberAblationVoluminaLabel->setText(QString::number(m_AblationZoneCentersProcessed.size()));
  mitk::RenderingManager::GetInstance()->Modified();
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();

  if (AblationUtils::CheckImageForNonAblatedTissue(m_SegmentationImage, m_ImageDimension) )
  {
    MITK_WARN << "There is still non ablated tumor tissue.";
  }
  std::chrono::milliseconds end = std::chrono::duration_cast< std::chrono::milliseconds >(
    std::chrono::system_clock::now().time_since_epoch());
  std::chrono::milliseconds diff = end - start;
  m_Stats.time = diff.count();
  this->CreateSpheresOfAblationVolumes();
  this->FillComboBoxAblationZones();
  this->CalculateAblationStatistics();
  this->SaveResults();
}

void QmitkAblationPlanningView::OnAblationRadiusChanged(double radius)
{
  if (m_AblationCalculationMade)
  {
    m_Controls.refreshCalculationsPushButton->setVisible(true);
  }
  m_AblationRadius = radius * (1 + ((double)m_Controls.tissueShrinkingSpinBox->value() / 100.0));
  MITK_DEBUG << "TissueShrinkingFactor increased ablation radius to: " << m_AblationRadius;
}

void QmitkAblationPlanningView::OnTissueShrinkingFactorChanged(int tissueShrinking)
{
  if (m_AblationCalculationMade)
  {
    m_Controls.refreshCalculationsPushButton->setVisible(true);
  }
  m_AblationRadius = m_Controls.ablationRadiusSpinBox->value() * (1 + ((double)tissueShrinking/100.0));
  MITK_DEBUG << "TissueShrinkingFactor increased ablation radius to: " << m_AblationRadius;
}

void QmitkAblationPlanningView::OnConfirmNewPositionClicked()
{
  if( m_Controls.ablationZonesComboBox->count() > 0 && m_SegmentationImage.IsNotNull())
  {
    int index = m_Controls.ablationZonesComboBox->currentIndex();

    AblationUtils::RemoveAblationVolume(m_AblationZoneCentersProcessed.at(index), m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);

    //Get the actual marked position of the crosshair:
    mitk::Point3D newPositionInWorldCoordinates = this->GetRenderWindowPart()->GetSelectedPosition();
    itk::Index<3> newPositionInIndexCoordinates;

    //Calculate the index coordinates of the new position:
    m_SegmentationImage->GetGeometry()->WorldToIndex(newPositionInWorldCoordinates,
      newPositionInIndexCoordinates);

    m_AblationZoneCentersProcessed.at(index) = newPositionInIndexCoordinates;
    m_AblationZoneCenters.at(index) = newPositionInIndexCoordinates;

    AblationUtils::CalculateAblationVolume(newPositionInIndexCoordinates, m_SegmentationImage, m_AblationRadius, m_ImageSpacing, m_ImageDimension);

    this->DeleteAllSpheres();
    this->CreateSpheresOfAblationVolumes();
    this->CalculateAblationStatistics();
  }
}

void QmitkAblationPlanningView::OnDeleteChosenAblationZoneClicked()
{
  if (m_Controls.ablationZonesComboBox->count() > 0 && m_SegmentationImage.IsNotNull())
  {
    int index = m_Controls.ablationZonesComboBox->currentIndex();
    AblationUtils::RemoveAblationVolume(m_AblationZoneCentersProcessed.at(index), m_SegmentationImage, m_AblationRadius, m_ImageDimension, m_ImageSpacing);

    this->DeleteAllSpheres();

    std::vector<itk::Index<3>>::iterator it = m_AblationZoneCentersProcessed.begin();
    m_AblationZoneCentersProcessed.erase(it + index);
    std::vector<itk::Index<3>>::iterator it2 = m_AblationZoneCenters.begin();
    m_AblationZoneCenters.erase(it2 + index);

    this->FillComboBoxAblationZones();
    this->CreateSpheresOfAblationVolumes();
    this->CalculateAblationStatistics();
  }
}

void QmitkAblationPlanningView::OnAddNewAblationZoneClicked()
{
  if( m_SegmentationImage.IsNotNull() )
  {
    if( m_TumorTissueSafetyMarginIndices.size() == 0 )
    {
      AblationUtils::FillVectorContainingIndicesOfTumorTissueSafetyMargin(m_SegmentationImage, m_ImageDimension, m_TumorTissueSafetyMarginIndices);
    }
    //Get the actual marked position of the crosshair:
    mitk::Point3D newPositionInWorldCoordinates = this->GetRenderWindowPart()->GetSelectedPosition();
    itk::Index<3> newPositionInIndexCoordinates;

    //Calculate the index coordinates of the new position:
    m_SegmentationImage->GetGeometry()->WorldToIndex(newPositionInWorldCoordinates,
      newPositionInIndexCoordinates);

    m_AblationZoneCentersProcessed.push_back(newPositionInIndexCoordinates);
    m_AblationZoneCenters.push_back(newPositionInIndexCoordinates);

    AblationUtils::CalculateAblationVolume(newPositionInIndexCoordinates, m_SegmentationImage, m_AblationRadius, m_ImageSpacing, m_ImageDimension);

    this->FillComboBoxAblationZones();
    this->DeleteAllSpheres();
    this->CreateSpheresOfAblationVolumes();
    this->CalculateAblationStatistics();
  }
}

void QmitkAblationPlanningView::OnCalculationModelChanged(bool)
{
  if (m_SegmentationImage.IsNotNull())
  {
    AblationUtils::ResetSegmentationImage(m_SegmentationImage, m_ImageDimension);
    this->DeleteAllSpheres();

    m_TumorTissueSafetyMarginIndices.clear();
    m_AblationZoneCenters.clear();
    m_TempAblationZoneCenters.clear();
    m_AblationZoneCentersProcessed.clear();
    m_TempAblationZoneCentersProcessed.clear();

    this->CalculateAblationStatistics();
    mitk::RenderingManager::GetInstance()->Modified();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
}

void QmitkAblationPlanningView::OnNumberOfRepetitionsChanged()
{
  if (m_AblationCalculationMade)
  {
    m_Controls.refreshCalculationsPushButton->setVisible(true);
  }
}

void QmitkAblationPlanningView::OnPercentageNonAblatedVolumeChanged()
{
  if (m_AblationCalculationMade)
  {
    m_Controls.refreshCalculationsPushButton->setVisible(true);
  }
}

void QmitkAblationPlanningView::SaveResults() {
  QString filename = m_Controls.outFolder->text() + m_Controls.expName->text();

  //Save CSV file with output
  QFile file(filename + ".csv");
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate))
  {
    MITK_WARN << "Cannot open file '" << filename.toStdString() << ".csv" << "' for writing.";
    return;
  }
  QTextStream outStream(&file);
  //write header
  outStream << "numberOfZones;tumorVolume[ml];safetyMarginVolume[ml];tumorAndSafetyMarginVolume[ml];totalAblationVolume[ml];ablationVolumeAblatedMoreThanOneTime[ml];factorOverlappingAblationZones;factorAblatedVolumeOutsideSafetyMargin;time[ms]\n";
  //write data
  outStream << m_Stats.numberOfZones << ";"
            << m_Stats.tumorVolume << ";"
            << m_Stats.safetyMarginVolume << ";"
            << m_Stats.tumorAndSafetyMarginVolume << ";"
            << m_Stats.totalAblationVolume << ";"
            << m_Stats.ablationVolumeAblatedMoreThanOneTime << ";"
            << m_Stats.factorOverlappingAblationZones << ";"
            << m_Stats.factorAblatedVolumeOutsideSafetyMargin << ";"
            << m_Stats.time << ";" << "\n";


  //Save MITK scene
  mitk::SceneIO::Pointer mySceneIO = mitk::SceneIO::New();
  QString filenameScene = filename + "_mitkScene.mitk";
  mitk::NodePredicateNot::Pointer isNotHelperObject =
    mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true)));
  mitk::DataStorage::SetOfObjects::ConstPointer nodesToBeSaved = this->GetDataStorage()->GetSubset(isNotHelperObject);
  mySceneIO->SaveScene(nodesToBeSaved, this->GetDataStorage(), filenameScene.toStdString().c_str());
}


void QmitkAblationPlanningView::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  OnTissueShrinkingFactorChanged(m_Controls.tissueShrinkingSpinBox->value());

  m_Controls.refreshCalculationsPushButton->setVisible(false);

  m_Controls.segmentationComboBox->SetDataStorage(GetDataStorage());
  m_Controls.segmentationComboBox->SetPredicate(m_IsASegmentationImagePredicate);

  // create signal/slot connections
  connect(m_Controls.segmentationComboBox, SIGNAL(OnSelectionChanged(const mitk::DataNode*)),
    this, SLOT(OnSegmentationComboBoxSelectionChanged(const mitk::DataNode*)));
  connect(m_Controls.ablationStartingPointPushButton, SIGNAL(clicked()),
    this, SLOT(OnAblationStartingPointPushButtonClicked()));
  connect(m_Controls.calculateAblationZonesPushButton, SIGNAL(clicked()),
    this, SLOT(OnCalculateAblationZonesPushButtonClicked()));
  connect(m_Controls.ablationRadiusSpinBox, SIGNAL(valueChanged(double)),
    this, SLOT(OnAblationRadiusChanged(double)));
  connect(m_Controls.tissueShrinkingSpinBox, SIGNAL(valueChanged(int)),
    this, SLOT(OnTissueShrinkingFactorChanged(int)));
  connect(m_Controls.calculateSafetyMarginPushButton, SIGNAL(clicked()),
    this, SLOT(OnCalculateSafetyMargin()));
  connect(m_Controls.confirmNewPositionPushButton, SIGNAL(clicked()),
    this, SLOT(OnConfirmNewPositionClicked()));
  connect(m_Controls.deleteChosenAblationZonePushButton, SIGNAL(clicked()),
    this, SLOT(OnDeleteChosenAblationZoneClicked()));
  connect(m_Controls.addNewAblationZonePushButton, SIGNAL(clicked()),
    this, SLOT(OnAddNewAblationZoneClicked()));
  connect(m_Controls.gridModelRadioButton, SIGNAL(toggled(bool)),
    this, SLOT(OnCalculationModelChanged(bool)));
  connect(m_Controls.randomDistributionRadioButton, SIGNAL(toggled(bool)),
    this, SLOT(OnCalculationModelChanged(bool)));
  connect(m_Controls.refreshCalculationsPushButton, SIGNAL(clicked()),
    this, SLOT(OnCalculateAblationZonesPushButtonClicked()));
  connect(m_Controls.repititionsCalculatingAblationZonesSpinBox, SIGNAL(valueChanged(int)),
    this, SLOT(OnNumberOfRepetitionsChanged()));
  connect(m_Controls.toleranceNonAblatedTumorSafetyMarginVolumeSpinBox, SIGNAL(valueChanged(int)),
    this, SLOT(OnPercentageNonAblatedVolumeChanged()));
  connect(m_Controls.mergeTumor, SIGNAL(clicked()), this, SLOT(OnMergeTumorSafetyMargin()));



  mitk::DataStorage::SetOfObjects::ConstPointer segmentationImages = GetDataStorage()->GetSubset(m_IsASegmentationImagePredicate);
  if (!segmentationImages->empty())
  {
    OnSelectionChanged(*segmentationImages->begin());
  }

  // set callback function for already existing nodes (segmentations)
  mitk::DataStorage::SetOfObjects::ConstPointer allImages = GetDataStorage()->GetSubset(m_IsOfTypeImagePredicate);
  for (mitk::DataStorage::SetOfObjects::const_iterator iter = allImages->begin(); iter != allImages->end(); ++iter)
  {
    mitk::DataNode* node = *iter;
    itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command = itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
    command->SetCallbackFunction(this, &QmitkAblationPlanningView::OnVisiblePropertyChanged);
    m_WorkingDataObserverTags.insert(std::pair<mitk::DataNode*, unsigned long>(node, node->GetProperty("visible")->AddObserver(itk::ModifiedEvent(), command)));

    itk::SimpleMemberCommand<QmitkAblationPlanningView>::Pointer command2 = itk::SimpleMemberCommand<QmitkAblationPlanningView>::New();
    command2->SetCallbackFunction(this, &QmitkAblationPlanningView::OnBinaryPropertyChanged);
    m_BinaryPropertyObserverTags.insert(std::pair<mitk::DataNode*, unsigned long>(node, node->GetProperty("binary")->AddObserver(itk::ModifiedEvent(), command2)));
  }

  /*itk::SimpleMemberCommand<QmitkSegmentationView>::Pointer command = itk::SimpleMemberCommand<QmitkSegmentationView>::New();
  command->SetCallbackFunction(this, &QmitkSegmentationView::RenderingManagerReinitialized);
  m_RenderingManagerObserverTag = mitk::RenderingManager::GetInstance()->AddObserver(mitk::RenderingManagerViewsInitializedEvent(), command);

  InitToolManagerSelection(m_Controls->patImageSelector->GetSelectedNode(), m_Controls->segImageSelector->GetSelectedNode());
  */
  /*m_RenderWindowPart = GetRenderWindowPart();
  if (m_RenderWindowPart)
  {
    RenderWindowPartActivated(m_RenderWindowPart);
  }*/
}
