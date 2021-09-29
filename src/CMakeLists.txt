
paraview_add_plugin(pvTipsyReader
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES TipsyReader
  MODULE_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")
