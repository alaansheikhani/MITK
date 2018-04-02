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

// mitk image
#include <mitkImage.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkLabelSetImage.h>
#include "mitkProperties.h"




const std::string QmitkAblationPlanningView::VIEW_ID = "org.mitk.views.ablationplanning";


//=====================Konstruktor/Destruktor===================================
QmitkAblationPlanningView::QmitkAblationPlanningView()
  :
  m_MouseCursorSet(false),
  m_DataSelectionChanged(false),
  m_AblationStartingPositionInWorldCoordinates(),
  m_AblationStartingPositionIndexCoordinates(),
  m_AblationStartingPositionValid(false)
{
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
  MITK_INFO << "OnSelectionChanged()";
  berry::IWorkbenchPart::Pointer nullPart;
  QList<mitk::DataNode::Pointer> nodes;
  nodes.push_back(node);
  this->OnSelectionChanged(nullPart, nodes);
}

void QmitkAblationPlanningView::OnSelectionChanged(berry::IWorkbenchPart::Pointer part, const QList<mitk::DataNode::Pointer>& nodes)
{
  MITK_INFO << "OnSelectionChanged()";
  if (nodes.size() != 0)
  {
    std::string markerName = "Position";
    unsigned int numberOfNodes = nodes.size();
    std::string nodeName = nodes.at(0)->GetName();
    if ((numberOfNodes == 1) && (nodeName.find(markerName) == 0))
    {
      //OnContourMarkerSelected(nodes.at(0));
      return;
    }
  }

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
        // set all nodes to invisible
        mitk::DataStorage::SetOfObjects::ConstPointer allImages = GetDataStorage()->GetSubset(m_IsOfTypeImagePredicate);
        for (mitk::DataStorage::SetOfObjects::const_iterator iter = allImages->begin(); iter != allImages->end(); ++iter)
        {
          (*iter)->SetVisibility(false);
        }

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
        /*else
        {
          // did not find a source / patient image, check all images and compare geometry
          mitk::DataStorage::SetOfObjects::ConstPointer possiblePatientImages = GetDataStorage()->GetSubset(m_IsAPatientImagePredicate);
          for (mitk::DataStorage::SetOfObjects::ConstIterator iter = possiblePatientImages->Begin(); iter != possiblePatientImages->End(); ++iter)
          {
            refNode = iter->Value();

            if (CheckForSameGeometry(selectedNode, iter->Value()))
            {
              refNode->SetVisibility(true);
              selectedNode->SetVisibility(true);
              SetToolManagerSelection(refNode, selectedNode);

              // doing this we can assure that the segmentation is always visible if the segmentation and the patient image are at the
              // same level in the data manager
              int layer(10);
              refNode->GetIntProperty("layer", layer);
              layer++;
              selectedNode->SetProperty("layer", mitk::IntProperty::New(layer));
              return;
            }
          }
          // did not find a source / patient image with the same geometry
          SetToolManagerSelection(nullptr, selectedNode);
        }*/
        mitk::RenderingManager::GetInstance()->InitializeViews(selectedNode->GetData()->GetTimeGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true);
      }
      else
      {
        MITK_WARN << "SelectedNode is no segmentation node";
      }
    }

}

//==============================================================================

void QmitkAblationPlanningView::SetFocus()
{
  m_Controls.segmentationComboBox->setFocus();
}

void QmitkAblationPlanningView::NodeRemoved(const mitk::DataNode * node)
{
  MITK_INFO << "NodeRemoved()";
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
  MITK_INFO << "NodeAdded()";
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

  //ApplyDisplayOptions(const_cast<mitk::DataNode*>(node));
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

void QmitkAblationPlanningView::OnSegmentationComboBoxSelectionChanged(const mitk::DataNode* node)
{
  if (node == nullptr)
  {
    //this->UpdateWarningLabel(tr("Select or create a segmentation"));
    //this->SetToolSelectionBoxesEnabled(false);
    return;
  }

  mitk::DataNode* refNode = m_Controls.segmentationComboBox->GetSelectedNode();

  //RenderingManagerReinitialized();
  //if (m_Controls->lblSegmentationWarnings->isVisible()) // "RenderingManagerReinitialized()" caused a warning. we do not need to go any further
  //  return;

  mitk::DataStorage::SetOfObjects::ConstPointer possibleParents = this->GetDataStorage()->GetSources(node, m_IsAPatientImagePredicate);

  if (possibleParents->Size() == 1)
  {
    mitk::DataNode* parentNode = possibleParents->ElementAt(0);

    if (parentNode != refNode)
    {
      //this->UpdateWarningLabel(tr("The selected segmentation does not match with the selected patient image!"));
      //this->SetToolSelectionBoxesEnabled(false);
      //this->SetToolManagerSelection(nullptr, node);
    }
    else
    {
      //this->UpdateWarningLabel("");
      //this->SetToolManagerSelection(refNode, node);
    }
  }
  else if (refNode && this->CheckForSameGeometry(node, refNode))
  {
    //this->UpdateWarningLabel("");
    //this->SetToolManagerSelection(refNode, node);
  }
  else if (!refNode || !this->CheckForSameGeometry(node, refNode))
  {
    //this->UpdateWarningLabel(tr("Please select or load the according patient image!"));
  }


  /*mitk::IRenderWindowPart* renderWindowPart = this->GetRenderWindowPart();
  if (!renderWindowPart || !node->IsVisible(renderWindowPart->GetQmitkRenderWindow("axial")->GetRenderer()))
  {
    this->UpdateWarningLabel(tr("The selected segmentation is currently not visible!"));
    this->SetToolSelectionBoxesEnabled(false);
  }*/
}

void QmitkAblationPlanningView::OnVisiblePropertyChanged()
{
  MITK_INFO << "OnVisiblePropertyChanged()";
  mitk::DataNode* selectedNode = m_Controls.segmentationComboBox->GetSelectedNode();
  if (!selectedNode)
  {
    return;
  }

  //mitk::IRenderWindowPart* renderWindowPart = this->GetRenderWindowPart();
  //bool selectedNodeIsVisible = renderWindowPart && selectedNode->IsVisible(renderWindowPart->GetQmitkRenderWindow("axial")->GetRenderer());

  /*if (!selectedNodeIsVisible)
  {
    this->SetToolSelectionBoxesEnabled(false);
    this->UpdateWarningLabel(tr("The selected segmentation is currently not visible!"));
  }
  else
  {
    this->SetToolSelectionBoxesEnabled(true);
    this->UpdateWarningLabel("");
  }*/
}

void QmitkAblationPlanningView::OnBinaryPropertyChanged()
{
  MITK_INFO << "OnBinaryPropertyChanged()";
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

void QmitkAblationPlanningView::OnAblationStatingPointPushButtonClicked()
{

  mitk::DataNode* selectedSegmentation = m_Controls.segmentationComboBox->GetSelectedNode();
  if( selectedSegmentation == nullptr )
  {
    m_AblationStartingPositionValid = false;
    return;
  }


  mitk::Image::Pointer segmentationImage = dynamic_cast<mitk::Image*>(selectedSegmentation->GetData());
  if (segmentationImage.IsNull())
  {
    MITK_WARN << "Failed to cast selected segmentation node to mitk::Image*";
    m_AblationStartingPositionValid = false;
    return;
  }

  //Get the actual marked position of the crosshair:
  m_AblationStartingPositionInWorldCoordinates = this->GetRenderWindowPart()->GetSelectedPosition();

  //Calculate the index coordinates of the starting position:
  segmentationImage->GetGeometry()->WorldToIndex( m_AblationStartingPositionInWorldCoordinates,
                                                  m_AblationStartingPositionIndexCoordinates);

  double pixelType = segmentationImage->GetPixelValueByIndex(m_AblationStartingPositionIndexCoordinates);
  MITK_INFO << "PixelType: " << pixelType;
  if (pixelType < 1.0)
  {
    m_AblationStartingPositionValid = false;
    m_Controls.ablationStartingPointLabel->setText("Position is not in the segmentation. Please choose a new starting position.");
    return;
  }

  m_AblationStartingPositionValid = true;
  double x = m_AblationStartingPositionInWorldCoordinates[0];
  double y = m_AblationStartingPositionInWorldCoordinates[1];
  double z = m_AblationStartingPositionInWorldCoordinates[2];
  QString text = QString("Set Ablation Startingposition to: %1 | %2 | %3").arg(x).arg(y).arg(z);
  m_Controls.ablationStartingPointLabel->setText(text);
  MITK_INFO << "Set Ablation Startingposition to: " << m_AblationStartingPositionInWorldCoordinates;
  MITK_INFO << "Startingposition in Index: " << m_AblationStartingPositionIndexCoordinates;
  MITK_INFO << "Spacing: " << segmentationImage->GetGeometry()->GetSpacing();
  //Get number of voxels in the three dimensions:
  MITK_INFO << "Dimension: " << segmentationImage->GetDimension(0) << " " << segmentationImage->GetDimension(1) << " " << segmentationImage->GetDimension(2);
}

void QmitkAblationPlanningView::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);

  m_Controls.segmentationComboBox->SetDataStorage(GetDataStorage());
  m_Controls.segmentationComboBox->SetPredicate(m_IsASegmentationImagePredicate);
  if (m_Controls.segmentationComboBox->GetSelectedNode().IsNotNull())
  {
    // TODO : update Text UpdateWarningLabel("");
  }

  // create signal/slot connections
  connect(m_Controls.segmentationComboBox, SIGNAL(OnSelectionChanged(const mitk::DataNode*)),
    this, SLOT(OnSegmentationComboBoxSelectionChanged(const mitk::DataNode*)));
  connect(m_Controls.ablationStartingPointPushButton, SIGNAL(clicked()),
    this, SLOT(OnAblationStatingPointPushButtonClicked()));


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
