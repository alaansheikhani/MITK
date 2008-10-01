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

#include "Poco/ClassLibrary.h"

#include <cherryIBundleActivator.h>
#include "cherryQtUIPlugin.h"

#include "cherryQtWorkbenchTweaklet.h"
#include "cherryQtWorkbenchPageTweaklet.h"
#include "internal/cherryQtWidgetsTweaklet.h"
#include "internal/defaultpresentation/cherryQtWorkbenchPresentationFactory.h"

#include <cherryIViewPart.h>
#include "cherryLogView.h"

POCO_BEGIN_MANIFEST(cherry::IBundleActivator)
  POCO_EXPORT_CLASS(cherry::QtUIPlugin)
POCO_END_MANIFEST

//****  Tweaklets  ****************

POCO_BEGIN_NAMED_MANIFEST(cherryWorkbenchTweaklet, cherry::WorkbenchTweaklet)
  POCO_EXPORT_CLASS(cherry::QtWorkbenchTweaklet)
POCO_END_MANIFEST

POCO_BEGIN_NAMED_MANIFEST(cherryWorkbenchPageTweaklet, cherry::WorkbenchPageTweaklet)
  POCO_EXPORT_CLASS(cherry::QtWorkbenchPageTweaklet)
POCO_END_MANIFEST

POCO_BEGIN_NAMED_MANIFEST(cherryGuiWidgetsTweaklet, cherry::GuiWidgetsTweaklet)
  POCO_EXPORT_CLASS(cherry::QtWidgetsTweaklet)
POCO_END_MANIFEST


//**** Views  **********************
POCO_BEGIN_NAMED_MANIFEST(cherryIViewPart, cherry::IViewPart)
  POCO_EXPORT_CLASS(cherry::LogView)
POCO_END_MANIFEST

//**** Presentation Factories *********
POCO_BEGIN_NAMED_MANIFEST(cherryIPresentationFactory, cherry::IPresentationFactory)
  POCO_EXPORT_CLASS(cherry::QtWorkbenchPresentationFactory)
POCO_END_MANIFEST
