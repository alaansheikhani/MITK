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


#ifndef CHERRYWORKBENCHTWEAKLET_H_
#define CHERRYWORKBENCHTWEAKLET_H_

#include "../internal/cherryTweaklets.h"

#include "../cherryShell.h"
#include "../dialogs/cherryIDialog.h"

namespace cherry
{

class WorkbenchWindow;

struct CHERRY_UI WorkbenchTweaklet : public Object
{
  cherryInterfaceMacro(WorkbenchTweaklet, cherry);

  static Tweaklets::TweakKey<WorkbenchTweaklet> KEY;

  static const std::string DIALOG_ID_SHOW_VIEW; // = "org.opencherry.ui.dialogs.showview";

  virtual int RunEventLoop() = 0;

  virtual bool IsRunning() = 0;

  /**
   * Creates the default contents and layout of the workbench window shell.
   * Returns the GUI dependent parent widget under which the workbench pages
   * create their contents.
   *
   * @param shell
   *            the shell
   */
  //virtual void* CreateDefaultWindowContents(Shell::Pointer shell) = 0;

  /**
   * Creates the client composite, under which workbench pages
   * create their controls.
   *
   */
  virtual void* CreatePageComposite(void* parent) = 0;

  virtual IDialog::Pointer CreateStandardDialog(const std::string& id) = 0;

};

}

#endif /* CHERRYWORKBENCHTWEAKLET_H_ */
