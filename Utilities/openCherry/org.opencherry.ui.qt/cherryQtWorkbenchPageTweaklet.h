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


#ifndef CHERRYQTWORKBENCHPAGETWEAKLET_H_
#define CHERRYQTWORKBENCHPAGETWEAKLET_H_

#include <tweaklets/cherryWorkbenchPageTweaklet.h>

#include "cherryUiQtDll.h"

namespace cherry
{

class CHERRY_UI_QT QtWorkbenchPageTweaklet : public WorkbenchPageTweaklet
{

public:

  cherryClassMacro(QtWorkbenchPageTweaklet);

  void* CreateClientComposite(void* pageControl);
  void* CreatePaneControl(void* parent);

  IViewPart::Pointer CreateErrorViewPart(const std::string& partName = "", const std::string& msg = "");
  IEditorPart::Pointer CreateErrorEditorPart(const std::string& partName = "", const std::string& msg = "");

};

}

#endif /* CHERRYQTWORKBENCHPAGETWEAKLET_H_ */
