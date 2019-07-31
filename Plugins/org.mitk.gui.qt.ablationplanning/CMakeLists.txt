project(org_mitk_gui_qt_ablationplanning)

mitk_create_plugin(
  EXPORT_DIRECTIVE ABLATIONPLANNING_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkMultilabel
)
