/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.s

============================================================================*/

#ifndef QmitknnUNetFolderParser_h_Included
#define QmitknnUNetFolderParser_h_Included

#include "QmitknnUNetToolGUI.h"
#include <QDirIterator>
#include <QString>
#include <vector>

struct FolderNode
{
  QString name;
  QString path; // parent
  std::vector<std::shared_ptr<FolderNode>> subFolders;
};

class MITKSEGMENTATIONUI_EXPORT QmitknnUNetFolderParser
{
public:

  /**
   * @brief Construct a new QmitknnUNetFolderParser object
   * Initializes root folder node object pointer calls
   * @param parentFolder 
   */
  QmitknnUNetFolderParser(const QString parentFolder)
  {
    m_RootNode = std::make_shared<FolderNode>();
    m_RootNode->path = parentFolder;
    m_RootNode->name = QString("nnUNet");
    m_RootNode->subFolders.clear();
    InitDirs(m_RootNode, 0);
  }

  /**
   * @brief Destroy the QmitknnUNetFolderParser object
   * 
   */
  ~QmitknnUNetFolderParser() = default; /*{ DeleteDirs(m_RootNode, LEVEL); }*/

  /**
   * @brief Returns the "Results Folder" string which is parent path of the root node.
   * 
   * @return QString 
   */
  QString getResultsFolder() { return m_RootNode->path; }

  /**
   * @brief Returns the Model Names from root node. Template function, 
   * type can be any of stl or Qt containers which supports push_back call.
   * 
   * @tparam T 
   * @return T 
   */
  template <typename T>
  T getModelNames()
  {
    auto models = GetSubFolderNamesFromNode<T>(m_RootNode);
    return models;
  }
  
  /**
   * @brief Returns the task names for a given model. Template function, 
   * type can be any of stl or Qt containers which supports push_back call.
   * 
   * @tparam T 
   * @return T 
   */
  template <typename T>
  T getTasksForModel(const QString &modelName)
  {
    std::shared_ptr<FolderNode> modelNode = GetSubNodeMatchingNameCrietria(modelName, m_RootNode);
    auto tasks = GetSubFolderNamesFromNode<T>(modelNode);
    return tasks;
  }

  /**
   * @brief Returns the models names for a given task. Template function, 
   * type can be any of stl or Qt containers which supports push_back call.
   * 
   * @tparam T 
   * @return T 
   */
  template <typename T>
  T getModelsForTask(const QString &taskName)
  {
    T modelsForTask;
    auto models = GetSubFolderNamesFromNode<T>(m_RootNode);
    foreach (QString model, models)
    {
      QStringList taskList = getTasksForModel<QStringList>(model);
      if (taskList.contains(taskName, Qt::CaseInsensitive))
      {
        modelsForTask << model;
      }
    }
    return modelsForTask;
  }

  /**
   * @brief Returns the trainer / planner names for a given task & model. Template function, 
   * type can be any of stl or Qt containers which supports push_back call.
   * 
   * @tparam T 
   * @return T 
   */
  template <typename T>
  T getTrainerPlannersForTask(const QString &taskName, const QString &modelName)
  {
    std::shared_ptr<FolderNode> modelNode = GetSubNodeMatchingNameCrietria(modelName, m_RootNode);
    std::shared_ptr<FolderNode> taskNode = GetSubNodeMatchingNameCrietria(taskName, modelNode);
    auto tps = GetSubFolderNamesFromNode<T>(taskNode);
    return tps;
  }

  /**
   * @brief Returns the Folds names for a given trainer,planner,task & model name. Template function, 
   * type can be any of stl or Qt containers which supports push_back call.
   * 
   * @tparam T 
   * @return T 
   */
  template <typename T>
  T getFoldsForTrainerPlanner(const QString &trainer,
                              const QString &planner,
                              const QString &taskName,
                              const QString &modelName)
  {
    std::shared_ptr<FolderNode> modelNode = GetSubNodeMatchingNameCrietria(modelName, m_RootNode);
    std::shared_ptr<FolderNode> taskNode = GetSubNodeMatchingNameCrietria(taskName, modelNode);
    QString trainerPlanner = trainer + QString("__") + planner;
    std::shared_ptr<FolderNode> tpNode = GetSubNodeMatchingNameCrietria(trainerPlanner, taskNode);
    auto folds = GetSubFolderNamesFromNode<T>(tpNode);
    return folds;
  }

private:
  const int m_LEVEL = 4;
  std::shared_ptr<FolderNode> m_RootNode;

  /**
   * @brief Iterates through the root node and returns the sub FolderNode object Matching Name Crietria 
   * 
   * @param queryName 
   * @param parentNode 
   * @return std::shared_ptr<FolderNode> 
   */
  std::shared_ptr<FolderNode> GetSubNodeMatchingNameCrietria(const QString &queryName,
                                                             std::shared_ptr<FolderNode> parentNode)
  {
    std::shared_ptr<FolderNode> retNode;
    std::vector<std::shared_ptr<FolderNode>> subNodes = parentNode->subFolders;
    for (std::shared_ptr<FolderNode> node : subNodes)
    {
      if (node->name == queryName)
      {
        retNode = node;
        break;
      }
    }
    return retNode;
  }

  /**
   * @brief Returns the sub folder names for a folder node object. Template function, 
   * type can be any of stl or Qt containers which supports push_back call.
   * 
   * @tparam T 
   * @return T 
   */
  template <typename T>
  T GetSubFolderNamesFromNode(const std::shared_ptr<FolderNode> parent)
  {
    T folders;
    std::vector<std::shared_ptr<FolderNode>> subNodes = parent->subFolders;
    for (std::shared_ptr<FolderNode> folder : subNodes)
    {
      folders.push_back(folder->name);
    }
    return folders;
  }

  /**
   * @brief Iterates through the sub folder hierarchy upto a level provided
   * and create a tree structure.
   * 
   * @param parent 
   * @param level 
   */
  void InitDirs(std::shared_ptr<FolderNode> parent, int level)
  {
    QString searchFolder = parent->path + QDir::separator() + parent->name;
    auto subFolders = FetchFoldersFromDir<QStringList>(searchFolder);
    level++;
    foreach (QString folder, subFolders)
    {
      std::shared_ptr<FolderNode> fp = std::make_shared<FolderNode>();
      fp->path = searchFolder;
      fp->name = folder;
      if (level < this->m_LEVEL)
      {
        InitDirs(fp, level);
      }
      parent->subFolders.push_back(fp);
    }
  }

  /**
   * @brief Iterates through the sub folder hierarchy upto a level provided
   * and clears the sub folder std::vector from each node.
   * 
   * @param parent 
   * @param level 
   */
  void DeleteDirs(std::shared_ptr<FolderNode> parent, int level)
  {
    level++;
    for (std::shared_ptr<FolderNode> subFolder : parent->subFolders)
    {
      if (level < m_LEVEL)
      {
        DeleteDirs(subFolder, level);
      }
      parent->subFolders.clear();
    }
  }

  /**
   * @brief Template function to fetch all folders inside a given path.
   * The type can be any of stl or Qt containers which supports push_back call.
   * 
   * @tparam T 
   * @param path 
   * @return T 
   */
  template <typename T>
  T FetchFoldersFromDir(const QString &path)
  {
    T folders;
    for (QDirIterator it(path, QDir::AllDirs, QDirIterator::NoIteratorFlags); it.hasNext();)
    {
      it.next();
      if (!it.fileName().startsWith('.'))
      {
        folders.push_back(it.fileName());
      }
    }
    return folders;
  }
};
#endif
