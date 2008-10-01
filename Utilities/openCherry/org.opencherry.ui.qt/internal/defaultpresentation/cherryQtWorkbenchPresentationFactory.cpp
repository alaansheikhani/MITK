/*=========================================================================

 Program:   openCherry Platform
 Language:  C++
 Date:      $Date$
 Version:   $Revision$

 Copyright (c) German Cancer Research Center, Division of Medical and
 Biological Informatics. All rights reserved.
 See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include "cherryQtWorkbenchPresentationFactory.h"

#include "cherryNativeTabFolder.h"
#include "cherryEmptyTabFolder.h"

#include "../util/cherryPresentablePartFolder.h"
#include "../util/cherryTabbedStackPresentation.h"

#include "../../cherryQtControlWidget.h"

namespace cherry
{

StackPresentation::Pointer QtWorkbenchPresentationFactory::CreateEditorPresentation(
    void* parent, IStackPresentationSite::Pointer site)
{
  NativeTabFolder* folder =
          new NativeTabFolder(static_cast<QWidget*>(parent));

  //    /*
  //     * Set the minimum characters to display, if the preference is something
  //     * other than the default. This is mainly intended for RCP applications
  //     * or for expert users (i.e., via the plug-in customization file).
  //     *
  //     * Bug 32789.
  //     */
  //    final IPreferenceStore store = PlatformUI.getPreferenceStore();
  //    if (store
  //        .contains(IWorkbenchPreferenceConstants.EDITOR_MINIMUM_CHARACTERS)) {
  //      final int minimumCharacters = store
  //          .getInt(IWorkbenchPreferenceConstants.EDITOR_MINIMUM_CHARACTERS);
  //      if (minimumCharacters >= 0) {
  //        folder.setMinimumCharacters(minimumCharacters);
  //      }
  //    }

  PresentablePartFolder* partFolder = new PresentablePartFolder(folder);

  StackPresentation::Pointer result =
          new TabbedStackPresentation(site, partFolder); //, new StandardEditorSystemMenu(site));

//  DefaultThemeListener themeListener =
//      new DefaultThemeListener(folder, result.getTheme());
//  result.getTheme().addListener(themeListener);
//
//  new DefaultMultiTabListener(result.getApiPreferences(), IWorkbenchPreferenceConstants.SHOW_MULTIPLE_EDITOR_TABS, folder);
//
//  new DefaultSimpleTabListener(result.getApiPreferences(), IWorkbenchPreferenceConstants.SHOW_TRADITIONAL_STYLE_TABS, folder);

  return result;
}

StackPresentation::Pointer QtWorkbenchPresentationFactory::CreateViewPresentation(
    void* parent, IStackPresentationSite::Pointer site)
{

  NativeTabFolder* folder =
          new NativeTabFolder(static_cast<QWidget*>(parent));

  //    final IPreferenceStore store = PlatformUI.getPreferenceStore();
  //    final int minimumCharacters = store
  //        .getInt(IWorkbenchPreferenceConstants.VIEW_MINIMUM_CHARACTERS);
  //    if (minimumCharacters >= 0) {
  //      folder.setMinimumCharacters(minimumCharacters);
  //    }

  PresentablePartFolder* partFolder = new PresentablePartFolder(folder);

  //folder->SetUnselectedCloseVisible(false);
  //folder->SetUnselectedImageVisible(true);

  StackPresentation::Pointer result =
          new TabbedStackPresentation(site, partFolder); //, new StandardViewSystemMenu(site));

//  DefaultThemeListener themeListener =
//      new DefaultThemeListener(folder, result.getTheme());
//  result.getTheme().addListener(themeListener);
//
//  new DefaultSimpleTabListener(result.getApiPreferences(), IWorkbenchPreferenceConstants.SHOW_TRADITIONAL_STYLE_TABS, folder);

  return result;
}

StackPresentation::Pointer QtWorkbenchPresentationFactory::CreateStandaloneViewPresentation(
    void* parent, IStackPresentationSite::Pointer site, bool showTitle)
{

  if (showTitle)
  {
    return this->CreateViewPresentation(parent, site);
  }
  EmptyTabFolder* folder = new EmptyTabFolder(static_cast<QWidget*>(parent), true);
  StackPresentation::Pointer presentation =
          new TabbedStackPresentation(site, folder); //, new StandardViewSystemMenu(site));

  return presentation;
}

std::string QtWorkbenchPresentationFactory::GetId()
{
  return "cherryQtWorkbenchPresentationFactory";
}

void* QtWorkbenchPresentationFactory::CreateSash(void* parent, int style)
{
  int orientation = style & (SASHORIENTATION_HORIZONTAL
      | SASHORIENTATION_VERTICAL);
  QWidget* sash = new QtControlWidget(static_cast<QWidget*>(parent));
  sash->setObjectName("Sash widget");
  if (orientation == SASHORIENTATION_HORIZONTAL)
    sash->setFixedHeight(this->GetSashSize(style));
  else
    sash->setFixedWidth(this->GetSashSize(style));

  return sash;
}

int QtWorkbenchPresentationFactory::GetSashSize(int style)
{
  return 3;
}

}
