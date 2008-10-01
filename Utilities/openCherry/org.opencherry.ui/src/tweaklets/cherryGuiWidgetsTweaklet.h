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


#ifndef CHERRYGUIWIDGETSTWEAKLET_H_
#define CHERRYGUIWIDGETSTWEAKLET_H_

#include <cherryMacros.h>

#include "../internal/cherryTweaklets.h"
#include "../guitk/cherryGuiTkISelectionListener.h"
#include "../guitk/cherryGuiTkIControlListener.h"
#include "../cherryRectangle.h"
#include "../cherryShell.h"

namespace cherry {

struct CHERRY_UI GuiWidgetsTweaklet : public Object
{
  cherryInterfaceMacro(GuiWidgetsTweaklet, cherry);

  static Tweaklets::TweakKey<GuiWidgetsTweaklet> KEY;

  virtual void AddSelectionListener(void* widget, GuiTk::ISelectionListener::Pointer listener) = 0;
  virtual void RemoveSelectionListener(void* widget, GuiTk::ISelectionListener::Pointer listener) = 0;

  /**
   * Adds the listener to the collection of listeners who will
   * be notified when the widget is moved or resized, by sending
   * it one of the messages defined in the <code>IControlListener</code>
   * interface.
   *
   * @param listener the listener which should be notified
   *
   * @see IControlListener
   * @see #RemoveControlListener
   */
  virtual void AddControlListener(void* widget, GuiTk::IControlListener::Pointer listener) = 0;

  /**
   * Removes the listener from the collection of listeners who will
   * be notified when the widget is moved or resized.
   *
   * @param listener the listener which should no longer be notified
   *
   * @see IControlListener
   * @see #AddControlListener
   */
  virtual void RemoveControlListener(void* widget, GuiTk::IControlListener::Pointer listener) = 0;

  virtual bool GetEnabled(void* widget) = 0;
  virtual void SetEnabled(void* widget, bool enabled) = 0;

  virtual void SetBounds(void* widget, const Rectangle& bounds) = 0;
  virtual Rectangle GetBounds(void* widget) = 0;

  virtual void SetVisible(void* widget, bool visible) = 0;
  virtual bool GetVisible(void* widget) = 0;
  virtual bool IsVisible(void* widget) = 0;

  virtual Rectangle GetClientArea(void* widget) = 0;

  virtual void* GetParent(void* widget) = 0;
  virtual bool SetParent(void* widget, void* parent) = 0;

  /**
   * Determines if one control is a child of another. Returns true iff the second
   * argument is a child of the first (or the same object).
   *
   * @param potentialParent
   * @param childToTest
   * @return
   */
  virtual bool IsChild(void* potentialParent, void* childToTest) = 0;

  /**
   * Returns the control which currently has keyboard focus,
   * or null if keyboard events are not currently going to
   * any of the controls built by the currently running
   * application.
   *
   * @return the control under the cursor
   */
  virtual void* GetFocusControl() = 0;

  virtual bool IsReparentable(void* widget) = 0;

  virtual void MoveAbove(void* widgetToMove, void* widget) = 0;
  virtual void MoveBelow(void* widgetToMove, void* widget) = 0;

  virtual void Dispose(void* widget) = 0;

  virtual Shell::Pointer CreateShell(Shell::Pointer parent, int style) = 0;
  virtual void DisposeShell(Shell::Pointer shell) = 0;

  virtual void* CreateComposite(void* parent) = 0;

  virtual std::vector<Shell::Pointer> GetShells() = 0;
  virtual Shell::Pointer GetShell(void* widget) = 0;
  virtual Shell::Pointer GetActiveShell() = 0;

  /**
   * Converts the given rectangle from display coordinates to the local coordinate system
   * of the given object
   *
   * @param coordinateSystem local coordinate system (widget) being converted to
   * @param toConvert rectangle to convert
   * @return a rectangle in control coordinates
   * @since 3.0
   */
  virtual Rectangle ToControl(void* coordinateSystem,
          const Rectangle& toConvert) = 0;

  /**
   * Converts the given rectangle from the local coordinate system of the given object
   * into display coordinates.
   *
   * @param coordinateSystem local coordinate system (widget) being converted from
   * @param toConvert rectangle to convert
   * @return a rectangle in display coordinates
   * @since 3.0
   */
  virtual Rectangle ToDisplay(void* coordinateSystem,
          const Rectangle& toConvert) = 0;
};

}

#endif /* CHERRYGUIWIDGETSTWEAKLET_H_ */
