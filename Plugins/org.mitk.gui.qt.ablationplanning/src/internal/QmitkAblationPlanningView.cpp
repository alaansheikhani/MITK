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
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>

#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <mitkSurface.h>

#include <cmath>

const std::string QmitkAblationPlanningView::VIEW_ID = "org.mitk.views.ablationplanning";
const static short ABLATION_VALUE = 2;
const static short TUMOR_NOT_YET_ABLATED = 1;
const static short NO_TUMOR_ISSUE = 0;
const static short SAFETY_MARGIN = 256;

//=====================Konstruktor/Destruktor===================================
QmitkAblationPlanningView::QmitkAblationPlanningView()
  :
  m_MouseCursorSet(false),
  m_DataSelectionChanged(false),
  m_AblationStartingPositionInWorldCoordinates(),
  m_AblationStartingPositionIndexCoordinates(),
  m_ManualAblationStartingPositionSet(false),
  m_AblationRadius(0.5)
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
      mitk::RenderingManager::GetInstance()->InitializeViews(selectedNode->GetData()->GetTimeGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true);
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

void QmitkAblationPlanningView::FillVectorContainingIndicesOfTumorTissueSafetyMargin()
{
  if( m_SegmentationImage.IsNotNull() )
  {
    MITK_INFO << "Detecting the tumor tissue and safety margin is in progress...";
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          if( imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED ||
              imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN )
          {
            m_TumorTissueSafetyMarginIndices.push_back(actualIndex);
          }
        }
      }
    }
  }
}

void QmitkAblationPlanningView::FindAblationStartingPosition()
{
  if( m_SegmentationImage.IsNotNull() )
  {
    MITK_INFO << "Finding a random ablation starting position...";
    int randomIndex1 = rand() % m_TumorTissueSafetyMarginIndices.size();
    int randomIndex2 = rand() % m_TumorTissueSafetyMarginIndices.size();
    int randomIndex3 = rand() % m_TumorTissueSafetyMarginIndices.size();
    int randomIndex4 = rand() % m_TumorTissueSafetyMarginIndices.size();
    int randomIndex5 = rand() % m_TumorTissueSafetyMarginIndices.size();
    MITK_INFO << "Gezogene Zufallszahl 1: " << randomIndex1;
    MITK_INFO << "Gezogene Zufallszahl 2: " << randomIndex2;
    MITK_INFO << "Gezogene Zufallszahl 3: " << randomIndex3;
    MITK_INFO << "Gezogene Zufallszahl 4: " << randomIndex4;
    MITK_INFO << "Gezogene Zufallszahl 5: " << randomIndex5;

    itk::Index<3> startingPosition1 = m_TumorTissueSafetyMarginIndices.at(randomIndex1);
    itk::Index<3> startingPosition2 = m_TumorTissueSafetyMarginIndices.at(randomIndex2);
    itk::Index<3> startingPosition3 = m_TumorTissueSafetyMarginIndices.at(randomIndex3);
    itk::Index<3> startingPosition4 = m_TumorTissueSafetyMarginIndices.at(randomIndex4);
    itk::Index<3> startingPosition5 = m_TumorTissueSafetyMarginIndices.at(randomIndex5);

    std::vector<itk::Index<3>> startingPositions;
    startingPositions.push_back(startingPosition1);
    startingPositions.push_back(startingPosition2);
    startingPositions.push_back(startingPosition3);
    startingPositions.push_back(startingPosition4);
    startingPositions.push_back(startingPosition5);

    double radius1 = this->CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition1);
    double radius2 = this->CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition2);
    double radius3 = this->CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition3);
    double radius4 = this->CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition4);
    double radius5 = this->CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(startingPosition5);
    std::vector<double> radiusVector;
    std::vector<double>::iterator result;
    radiusVector.push_back(radius1);
    radiusVector.push_back(radius2);
    radiusVector.push_back(radius3);
    radiusVector.push_back(radius4);
    radiusVector.push_back(radius5);
    result = std::max_element(radiusVector.begin(), radiusVector.end());
    int index = std::distance(radiusVector.begin(), result);

    m_AblationStartingPositionIndexCoordinates = startingPositions.at(index);

    //Calculate the index coordinates of the starting position:
    m_SegmentationImage->GetGeometry()->IndexToWorld(m_AblationStartingPositionIndexCoordinates,
                                                     m_AblationStartingPositionInWorldCoordinates);

    double x = m_AblationStartingPositionInWorldCoordinates[0];
    double y = m_AblationStartingPositionInWorldCoordinates[1];
    double z = m_AblationStartingPositionInWorldCoordinates[2];
    QString text = QString("Set Ablation Startingposition to: %1 | %2 | %3").arg(x).arg(y).arg(z);
    m_Controls.ablationStartingPointLabel->setText(text);
    MITK_INFO << "Set Ablation Startingposition to: " << m_AblationStartingPositionInWorldCoordinates;
    MITK_INFO << "Startingposition in Index: " << m_AblationStartingPositionIndexCoordinates;
    MITK_INFO << "Spacing: " << m_SegmentationImage->GetGeometry()->GetSpacing();
    //Get number of voxels in the three dimensions:
    MITK_INFO << "Dimension: " << m_ImageDimension[0] << " " << m_ImageDimension[1] << " " << m_ImageDimension[2];
  }
}

double QmitkAblationPlanningView::CalculateScalarDistance(itk::Index<3> &point1, itk::Index<3> &point2)
{

  double x = (point1[0] - point2[0]) * m_ImageSpacing[0];
  double y = (point1[1] - point2[1]) * m_ImageSpacing[1];
  double z = (point1[2] - point2[2]) * m_ImageSpacing[2];

  return sqrt(x*x + y*y + z*z);
}

void QmitkAblationPlanningView::CalculateAblationVolume(itk::Index<3>& center)
{
  MITK_INFO << "Calculate ablation volume for index: " << center;
  if( m_SegmentationImage.IsNotNull() )
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for( actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1 )
    {
      for(actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for( actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          if( m_AblationRadius >= this->CalculateScalarDistance(center, actualIndex))
          {
            unsigned short pixelValue =
              imagePixelWriter.GetPixelByIndex(actualIndex)
                + ABLATION_VALUE;

            imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
          }
        }
      }
    }
    m_AblationZoneCenters.push_back(center);
  }
}

bool QmitkAblationPlanningView::CheckVolumeForNonAblatedTissue(itk::Index<3> &centerOfVolume)
{
  if(m_SegmentationImage.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(m_AblationRadius / m_ImageSpacing[0]);
    unsigned int pixelDirectionY = floor(m_AblationRadius / m_ImageSpacing[1]);
    unsigned int pixelDirectionZ = floor(m_AblationRadius / m_ImageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    this->CalculateUpperLowerXYZ(upperX, lowerX, upperY, lowerY, upperZ, lowerZ,
      pixelDirectionX, pixelDirectionY, pixelDirectionZ, centerOfVolume);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if (m_AblationRadius >= this->CalculateScalarDistance(centerOfVolume, actualIndex))
          {
            if( imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED ||
                imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN )
            {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

bool QmitkAblationPlanningView::CheckIfVolumeOfGivenRadiusIsTotallyInsideTumorTissueAndSafetyMargin(
                                double &radius, itk::Index<3>& centerOfVolume)
{
  if (m_SegmentationImage.IsNotNull())
  {
    unsigned int pixelDirectionX = floor(radius / m_ImageSpacing[0]);
    unsigned int pixelDirectionY = floor(radius / m_ImageSpacing[1]);
    unsigned int pixelDirectionZ = floor(radius / m_ImageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    this->CalculateUpperLowerXYZ(upperX, lowerX, upperY, lowerY, upperZ, lowerZ,
      pixelDirectionX, pixelDirectionY, pixelDirectionZ, centerOfVolume);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if( radius >= this->CalculateScalarDistance(centerOfVolume, actualIndex) )
          {
            if( imagePixelWriter.GetPixelByIndex(actualIndex) != TUMOR_NOT_YET_ABLATED &&
                imagePixelWriter.GetPixelByIndex(actualIndex) != SAFETY_MARGIN)
            {
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}

double QmitkAblationPlanningView::
       CalculateMaxRadiusOfVolumeInsideTumorForGivenPoint(itk::Index<3>& point)
{
  double radius = 1.0;
  while( this->CheckIfVolumeOfGivenRadiusIsTotallyInsideTumorTissueAndSafetyMargin(
                                                                      radius, point))
  {
    radius += 1;
  }

  MITK_INFO << "Calculated max radius for given point: " << radius;
  return radius;
}

bool QmitkAblationPlanningView::CheckImageForNonAblatedTissue()
{
  if(m_SegmentationImage.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          if( imagePixelWriter.GetPixelByIndex(actualIndex) == TUMOR_NOT_YET_ABLATED ||
              imagePixelWriter.GetPixelByIndex(actualIndex) == SAFETY_MARGIN)
          {
            return true;
          }
        }
      }
    }
  }
  return false;
}

void QmitkAblationPlanningView::ProcessDirectNeighbourAblationZones(itk::Index<3>& center)
{
  MITK_INFO << "Process direct neighbour ablation zones for index: " << center;
  std::vector<itk::Index<3>> indices =
    this->CalculateIndicesOfDirectNeighbourAblationZones(center);

  for( std::vector<itk::Index<3>>::iterator it = indices.begin();
       it != indices.end(); ++it )
  {
    if( CheckVolumeForNonAblatedTissue(*it) )
    {
      this->CalculateAblationVolume(*it);
    }
  }

  //Now, all 12 direct neighbour ablation zones are processed. So add the
  // index of the given center to the processed ablation centers:
  m_AblationZoneCentersProcessed.push_back(center);
}

void QmitkAblationPlanningView::CalculateUpperLowerXYZ( unsigned int &upperX,
                                                        unsigned int &lowerX,
                                                        unsigned int &upperY,
                                                        unsigned int &lowerY,
                                                        unsigned int &upperZ,
                                                        unsigned int &lowerZ,
                                                        unsigned int &pixelDirectionX,
                                                        unsigned int &pixelDirectionY,
                                                        unsigned int &pixelDirectionZ,
                                                        itk::Index<3> &center )
{
  //Calculate upperX --> means vector in direction [1,0,0]:
  if (center[0] + pixelDirectionX >= m_ImageDimension[0])
  {
    upperX = m_ImageDimension[0] - 1;
  }
  else
  {
    upperX = center[0] + pixelDirectionX;
  }
  //Calculate lowerX --> means vector in direction [-1,0,0]:
  if (center[0] - pixelDirectionX < 0)
  {
    lowerX = 0;
  }
  else
  {
    lowerX = center[0] - pixelDirectionX;
  }

  //Calculate upperY --> means vector in direction [0,1,0]:
  if (center[1] + pixelDirectionY >= m_ImageDimension[1])
  {
    upperY = m_ImageDimension[1] - 1;
  }
  else
  {
    upperY = center[1] + pixelDirectionY;
  }
  //Calculate lowerY --> means vector in direction [0,-1,0]:
  if (center[1] - pixelDirectionY < 0)
  {
    lowerY = 0;
  }
  else
  {
    lowerY = center[1] - pixelDirectionY;
  }

  //Calculate upperZ --> means vector in direction [0,0,1]:
  if (center[2] + pixelDirectionZ >= m_ImageDimension[2])
  {
    upperZ = m_ImageDimension[2] - 1;
  }
  else
  {
    upperZ = center[2] + pixelDirectionZ;
  }

  //Calculate lowerZ --> means vector in direction [0,0,-1]:
  if (center[2] - pixelDirectionZ < 0)
  {
    lowerZ = 0;
  }
  else
  {
    lowerZ = center[2] - pixelDirectionZ;
  }
}

std::vector<itk::Index<3>>
  QmitkAblationPlanningView::
  CalculateIndicesOfDirectNeighbourAblationZones(itk::Index<3>& center)
{
  MITK_INFO << "Calculate indices of direct neighbour ablation zones...";
  std::vector<itk::Index<3>> directNeighbourAblationZones;
  unsigned int pixelDirectionX = floor(m_AblationRadius / m_ImageSpacing[0]);
  unsigned int pixelDirectionY = floor(m_AblationRadius / m_ImageSpacing[1]);
  unsigned int pixelDirectionZ = floor(m_AblationRadius / m_ImageSpacing[2]);

  unsigned int upperX;
  unsigned int lowerX;
  unsigned int upperY;
  unsigned int lowerY;
  unsigned int upperZ;
  unsigned int lowerZ;

  itk::Index<3> newIndex;

  this->CalculateUpperLowerXYZ( upperX, lowerX, upperY, lowerY, upperZ, lowerZ,
                                pixelDirectionX, pixelDirectionY, pixelDirectionZ, center);

  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [1, 1, 0]:
  newIndex[0] = upperX;
  newIndex[1] = upperY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [1, 1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [-1, 1, 0]:
  newIndex[0] = lowerX;
  newIndex[1] = upperY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [-1, 1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [1, -1, 0]:
  newIndex[0] = upperX;
  newIndex[1] = lowerY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [1, -1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [-1, -1, 0]:
  newIndex[0] = lowerX;
  newIndex[1] = lowerY;
  newIndex[2] = center[2];

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [-1, -1, 0] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [1, 0, 1]:
  newIndex[0] = upperX;
  newIndex[1] = center[1];
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [1, 0, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [1, 0, -1]:
  newIndex[0] = upperX;
  newIndex[1] = center[1];
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [1, 0, -1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [-1, 0, 1]:
  newIndex[0] = lowerX;
  newIndex[1] = center[1];
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [-1, 0, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [-1, 0, -1]:
  newIndex[0] = lowerX;
  newIndex[1] = center[1];
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [-1, 0, -1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [0, 1, 1]:
  newIndex[0] = center[0];
  newIndex[1] = upperY;
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [0, 1, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [0, 1, -1]:
  newIndex[0] = center[0];
  newIndex[1] = upperY;
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [0, 1, -1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [0, -1, 1]:
  newIndex[0] = center[0];
  newIndex[1] = lowerY;
  newIndex[2] = upperZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [0, -1, 1] --> " << newIndex;
  //--------------------------------------------------------------------------------------
  //Calculate position in vector direction [0, -1, -1]:
  newIndex[0] = center[0];
  newIndex[1] = lowerY;
  newIndex[2] = lowerZ;

  directNeighbourAblationZones.push_back(newIndex);
  MITK_INFO << "Index for [0, -1, -1] --> " << newIndex;

  return directNeighbourAblationZones;
}

bool QmitkAblationPlanningView::IsAblationZoneAlreadyProcessed(itk::Index<3>& center)
{
  for( std::vector<itk::Index<3>>::iterator it = m_AblationZoneCentersProcessed.begin();
       it != m_AblationZoneCentersProcessed.end(); ++it )
  {
    if (center == (*it))
    {
      return true;
    }
  }
  return false;
}

void QmitkAblationPlanningView::DetectNotNeededAblationVolume()
{
  std::vector<int> indicesRemoved;
  for (int index = 0; index < m_AblationZoneCentersProcessed.size(); ++index)
  {
    if (!this->CheckIfAblationVolumeIsNeeded(m_AblationZoneCentersProcessed.at(index)))
    {
      this->RemoveAblationVolume(m_AblationZoneCentersProcessed.at(index));
      indicesRemoved.push_back(index);
    }
  }
  for( int index = indicesRemoved.size() - 1; index >= 0; --index )
  {
    std::vector<itk::Index<3>>::iterator it = m_AblationZoneCentersProcessed.begin();
    m_AblationZoneCentersProcessed.erase(it + indicesRemoved.at(index));
    std::vector<itk::Index<3>>::iterator it2 = m_AblationZoneCenters.begin();
    m_AblationZoneCenters.erase(it2 + indicesRemoved.at(index));
    MITK_INFO << "Removed Ablation zone at index position: " << indicesRemoved.at(index);
  }
}

bool QmitkAblationPlanningView::CheckIfAblationVolumeIsNeeded(itk::Index<3>& center)
{
  if( m_SegmentationImage.IsNotNull() )
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          if (m_AblationRadius >= this->CalculateScalarDistance(center, actualIndex))
          {
            if( imagePixelWriter.GetPixelByIndex(actualIndex) - ABLATION_VALUE == TUMOR_NOT_YET_ABLATED
                || imagePixelWriter.GetPixelByIndex(actualIndex) - ABLATION_VALUE == SAFETY_MARGIN )
            {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

void QmitkAblationPlanningView::RemoveAblationVolume(itk::Index<3>& center)
{
  if(m_SegmentationImage.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          if (m_AblationRadius >= this->CalculateScalarDistance(center, actualIndex))
          {
            unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
            pixelValue -= ABLATION_VALUE;
            imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
          }
        }
      }
    }
  }
}

void QmitkAblationPlanningView::CreateSpheresOfAblationVolumes()
{
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
    vtkSphere->SetCenter(centerInWorldCoordinates[0], centerInWorldCoordinates[1], centerInWorldCoordinates[2]);
    vtkSphere->Update();

    appendPolyData->AddInputData(vtkSphere->GetOutput());

    mySphere->SetVtkPolyData(vtkSphere->GetOutput());

    m_DataNode->SetData(mySphere);
    QString name;
    if (index == 0)
    {
      name = QString("Start-Ablationsvolumen");
    }
    else
    {
      name = QString("Kugel_%1").arg(index + 1);
    }

    m_DataNode->SetName(name.toStdString());
    this->GetDataStorage()->Add(m_DataNode);
  }

  this->GetDataStorage()->Modified();
  this->RequestRenderWindowUpdate();
}

void QmitkAblationPlanningView::DeleteAllSpheres()
{
  for( int index = m_AblationZoneCentersProcessed.size(); index > 0; --index )
  {
    QString name;
    if( index > 1 )
    {
      name = QString("Kugel_%1").arg(index);
    }
    else
    {
      name = QString("Start-Ablationsvolumen");
    }

    mitk::DataNode::Pointer dataNode = this->GetDataStorage()->GetNamedNode(name.toStdString());
    if( dataNode.IsNotNull() )
    {
      this->GetDataStorage()->Remove(dataNode);
    }
  }
  this->GetDataStorage()->Modified();
  this->RequestRenderWindowUpdate();
}


void QmitkAblationPlanningView::ResetSegmentationImage()
{
  if( m_SegmentationImage.IsNotNull() )
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
            unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
            pixelValue &= (SAFETY_MARGIN + TUMOR_NOT_YET_ABLATED);
            imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
        }
      }
    }
  }
}

void QmitkAblationPlanningView::ResetSafetyMargin()
{
  if (m_SegmentationImage.IsNotNull())
  {
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          unsigned short pixelValue = imagePixelWriter.GetPixelByIndex(actualIndex);
          pixelValue &= (SAFETY_MARGIN - 1);
          imagePixelWriter.SetPixelByIndex(actualIndex, pixelValue);
        }
      }
    }
  }
}

bool QmitkAblationPlanningView::CheckAllVonNeumannNeighbourPixelsAreTumorTissue(itk::Index<3> &pixel)
{
  if (m_SegmentationImage.IsNotNull())
  {
    //MITK_INFO << "CheckAllVonNeumannNeighbourPixels... " << pixel[0] << " " << pixel[1] << " " << pixel[2];
    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    if (imagePixelWriter.GetPixelByIndex(pixel) != TUMOR_NOT_YET_ABLATED)
    {
      //If true --> actual pixel is no tumor tissue, so skip this pixel:
      return true;
    }

    unsigned int xLower = pixel[0];
    if( xLower != 0 )
      --xLower;

    unsigned int xUpper = pixel[0];
    if( xUpper != m_ImageDimension[0] - 1 )
      ++xUpper;

    unsigned int yLower = pixel[1];
    if (yLower != 0)
      --yLower;

    unsigned int yUpper = pixel[1];
    if (yUpper != m_ImageDimension[1] - 1)
      ++yUpper;

    unsigned int zLower = pixel[2];
    if (zLower != 0)
      --zLower;

    unsigned int zUpper = pixel[2];
    if (zUpper != m_ImageDimension[2] - 1)
      ++zUpper;

    itk::Index<3> pixelXLower = pixel;
    pixelXLower[0] = xLower;
    itk::Index<3> pixelXUpper = pixel;
    pixelXUpper[0] = xUpper;

    itk::Index<3> pixelYLower = pixel;
    pixelYLower[1] = yLower;
    itk::Index<3> pixelYUpper = pixel;
    pixelYUpper[1] = yUpper;

    itk::Index<3> pixelZLower = pixel;
    pixelZLower[2] = zLower;
    itk::Index<3> pixelZUpper = pixel;
    pixelZUpper[2] = zUpper;

    if( imagePixelWriter.GetPixelByIndex(pixelXLower) != TUMOR_NOT_YET_ABLATED )
    {
      return false;
    }
    if (imagePixelWriter.GetPixelByIndex(pixelYLower) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if (imagePixelWriter.GetPixelByIndex(pixelZLower) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if (imagePixelWriter.GetPixelByIndex(pixelXUpper) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if (imagePixelWriter.GetPixelByIndex(pixelYUpper) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
    if (imagePixelWriter.GetPixelByIndex(pixelZUpper) != TUMOR_NOT_YET_ABLATED)
    {
      return false;
    }
  }
  return true;
}

void QmitkAblationPlanningView::CreateSafetyMarginInfluenceAreaOfPixel(itk::Index<3>& pixel)
{
  if( m_SegmentationImage.IsNotNull() )
  {
    double margin = m_Controls.safetyMarginSpinBox->value();
    unsigned int pixelDirectionX = floor(margin / m_ImageSpacing[0]);
    unsigned int pixelDirectionY = floor(margin / m_ImageSpacing[1]);
    unsigned int pixelDirectionZ = floor(margin / m_ImageSpacing[2]);

    unsigned int upperX;
    unsigned int lowerX;
    unsigned int upperY;
    unsigned int lowerY;
    unsigned int upperZ;
    unsigned int lowerZ;

    this->CalculateUpperLowerXYZ(upperX, lowerX, upperY, lowerY, upperZ, lowerZ,
       pixelDirectionX, pixelDirectionY, pixelDirectionZ, pixel);

    mitk::ImagePixelWriteAccessor<unsigned short, 3> imagePixelWriter(m_SegmentationImage);
    itk::Index<3> actualIndex;
    for (actualIndex[2] = lowerZ; actualIndex[2] <= upperZ; actualIndex[2] += 1)
    {
      for (actualIndex[1] = lowerY; actualIndex[1] <= upperY; actualIndex[1] += 1)
      {
        for (actualIndex[0] = lowerX; actualIndex[0] <= upperX; actualIndex[0] += 1)
        {
          if( margin >= this->CalculateScalarDistance(pixel, actualIndex) )
          {
            if( imagePixelWriter.GetPixelByIndex(actualIndex) != TUMOR_NOT_YET_ABLATED &&
                imagePixelWriter.GetPixelByIndex(actualIndex) != SAFETY_MARGIN )
            {
              imagePixelWriter.SetPixelByIndex(actualIndex, SAFETY_MARGIN );
            }
          }
        }
      }
    }
  }
}

void QmitkAblationPlanningView::OnSegmentationComboBoxSelectionChanged(const mitk::DataNode* node)
{
  MITK_INFO << "OnSegmentationComboBoxSelectionChanged()";
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
  MITK_INFO << "OnVisiblePropertyChanged()";

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

void QmitkAblationPlanningView::OnCalculateSafetyMargin()
{
  if(m_SegmentationImage.IsNotNull() && m_Controls.safetyMarginSpinBox->value() > 0.0 )
  {
    MITK_INFO << "Calculating safety margin is in progress...";
    //Reset the old calculated safety margin if the user calculated a safety margin before:
    this->ResetSafetyMargin();
    itk::Index<3> actualIndex;
    for (actualIndex[2] = 0; actualIndex[2] < m_ImageDimension[2]; actualIndex[2] += 1)
    {
      for (actualIndex[1] = 0; actualIndex[1] < m_ImageDimension[1]; actualIndex[1] += 1)
      {
        for (actualIndex[0] = 0; actualIndex[0] < m_ImageDimension[0]; actualIndex[0] += 1)
        {
          // Check, if von-Neumann neighbour elements are tumor tissue or not:
          if( !this->CheckAllVonNeumannNeighbourPixelsAreTumorTissue(actualIndex) )
          {
            this->CreateSafetyMarginInfluenceAreaOfPixel(actualIndex);
          }
        }
      }
    }
    MITK_INFO << "Finished calculating safety margin.";
    mitk::RenderingManager::GetInstance()->Modified();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
  else if(m_Controls.safetyMarginSpinBox->value() > 0.0)
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

  MITK_INFO << "PixelType: " << pixelType;
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
  MITK_INFO << "Set Ablation Startingposition to: " << m_AblationStartingPositionInWorldCoordinates;
  MITK_INFO << "Startingposition in Index: " << m_AblationStartingPositionIndexCoordinates;
  MITK_INFO << "Spacing: " << m_SegmentationImage->GetGeometry()->GetSpacing();
  //Get number of voxels in the three dimensions:
  MITK_INFO << "Dimension: " << m_ImageDimension[0] << " " << m_ImageDimension[1] << " " << m_ImageDimension[2];

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

  //Reset earlier calculations:
  this->DeleteAllSpheres();
  this->ResetSegmentationImage();
  m_TumorTissueSafetyMarginIndices.clear();
  m_AblationZoneCenters.clear();
  m_AblationZoneCentersProcessed.clear();

  this->FillVectorContainingIndicesOfTumorTissueSafetyMargin();

  if (!m_ManualAblationStartingPositionSet)
  {
    this->FindAblationStartingPosition();
  }

  this->CalculateAblationVolume(m_AblationStartingPositionIndexCoordinates);
  this->ProcessDirectNeighbourAblationZones(m_AblationStartingPositionIndexCoordinates);
  while( m_AblationZoneCenters.size() != m_AblationZoneCentersProcessed.size() )
  {
    MITK_INFO << "Size1: " << m_AblationZoneCenters.size() << " Size2: " << m_AblationZoneCentersProcessed.size();
    for( int index = 0; index < m_AblationZoneCenters.size(); ++index )
    {
      if (!this->IsAblationZoneAlreadyProcessed(m_AblationZoneCenters.at(index)))
      {
        this->ProcessDirectNeighbourAblationZones(m_AblationZoneCenters.at(index));
        break;
      }
    }

  }
  MITK_INFO << "Finished calculating ablation zones!";
  MITK_INFO << "Total number of ablation zones: " << m_AblationZoneCentersProcessed.size();

  this->DetectNotNeededAblationVolume();

  m_Controls.numberAblationVoluminaLabel->setText(QString::number(m_AblationZoneCentersProcessed.size()));
  mitk::RenderingManager::GetInstance()->Modified();
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();

  if (this->CheckImageForNonAblatedTissue())
  {
    MITK_WARN << "There is still non ablated tumor tissue.";
  }

  this->CreateSpheresOfAblationVolumes();

}

void QmitkAblationPlanningView::OnAblationRadiusChanged(double radius)
{
  m_AblationRadius = radius;
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
    this, SLOT(OnAblationStartingPointPushButtonClicked()));
  connect(m_Controls.calculateAblationZonesPushButton, SIGNAL(clicked()),
    this, SLOT(OnCalculateAblationZonesPushButtonClicked()));
  connect(m_Controls.ablationRadiusSpinBox, SIGNAL(valueChanged(double)),
    this, SLOT(OnAblationRadiusChanged(double)));
  connect(m_Controls.calculateSafetyMarginPushButton, SIGNAL(clicked()),
    this, SLOT(OnCalculateSafetyMargin()));

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
