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

#ifndef CHERRYWORKBENCHCONFIGURER_H_
#define CHERRYWORKBENCHCONFIGURER_H_

#include "../application/cherryIWorkbenchConfigurer.h"

namespace cherry
{

/**
 * Internal class providing special access for configuring the workbench.
 * <p>
 * Note that these objects are only available to the main application
 * (the plug-in that creates and owns the workbench).
 * </p>
 * <p>
 * This class is not intended to be instantiated or subclassed by clients.
 * </p>
 *
 * @since 3.0
 */
class WorkbenchConfigurer : public IWorkbenchConfigurer {

public:

  cherryClassMacro(WorkbenchConfigurer);

private:

    /**
     * Table to hold arbitrary key-data settings (key type: <code>String</code>,
     * value type: <code>Object</code>).
     * @see #setData
     */
    //Map extraData = new HashMap();

    /**
     * Indicates whether workbench state should be saved on close and
     * restored on subsequent open.
     */
    bool saveAndRestore;

    /**
     * Indicates whether the workbench is being force to close. During
     * an emergency close, no interaction with the user should be done.
     */
    bool isEmergencyClosing;

    /**
     * Indicates the behaviour when the last window is closed.
     * If <code>true</code>, the workbench will exit (saving the last window's state,
     * if configured to do so).
     * If <code>false</code> the window will be closed, leaving the workbench running.
     *
     * @since 3.1
     */
  bool exitOnLastWindowClose;

public:

    /**
     * Creates a new workbench configurer.
     * <p>
     * This method is declared package-private. Clients are passed an instance
     * only via {@link WorkbenchAdvisor#initialize WorkbenchAdvisor.initialize}
     * </p>
     */
    WorkbenchConfigurer();

    /* (non-javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#getWorkbench
     */
    IWorkbench* GetWorkbench();


    /* (non-javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#declareImage
     */
//    void declareImage(String symbolicName, ImageDescriptor descriptor,
//            boolean shared) {
//        if (symbolicName == null || descriptor == null) {
//            throw new IllegalArgumentException();
//        }
//        WorkbenchImages.declareImage(symbolicName, descriptor, shared);
//    }

    /* (non-javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#getWindowConfigurer
     */
    IWorkbenchWindowConfigurer::Pointer GetWindowConfigurer(IWorkbenchWindow::Pointer window);

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#getSaveAndRestore()
     */
    bool GetSaveAndRestore();

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#setSaveAndRestore(boolean)
     */
    void SetSaveAndRestore(bool enabled);

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#getData
     */
//    Object getData(String key) {
//        if (key == null) {
//            throw new IllegalArgumentException();
//        }
//        return extraData.get(key);
//    }

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#setData(String, Object)
     */
//    void setData(String key, Object data) {
//        if (key == null) {
//            throw new IllegalArgumentException();
//        }
//        if (data != null) {
//            extraData.put(key, data);
//        } else {
//            extraData.remove(key);
//        }
//    }

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#emergencyClose()
     */
    void EmergencyClose();

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#emergencyClosing()
     */
    bool EmergencyClosing();

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#restoreState()
     */
    bool RestoreState();

    /* (non-Javadoc)
     * @see org.opencherry.ui.application.IWorkbenchConfigurer#openFirstTimeWindow()
     */
    void OpenFirstTimeWindow();

  /* (non-Javadoc)
   * @see org.opencherry.ui.application.IWorkbenchConfigurer#restoreWorkbenchWindow(org.opencherry.ui.IMemento)
   */
  IWorkbenchWindowConfigurer::Pointer RestoreWorkbenchWindow(IMemento::Pointer memento);

  /* (non-Javadoc)
   * @see org.opencherry.ui.application.IWorkbenchConfigurer#getExitOnLastWindowClose()
   */
  bool GetExitOnLastWindowClose();

  /* (non-Javadoc)
   * @see org.opencherry.ui.application.IWorkbenchConfigurer#setExitOnLastWindowClose(boolean)
   */
  void SetExitOnLastWindowClose(bool enabled);
};

}

#endif /*CHERRYWORKBENCHCONFIGURER_H_*/
