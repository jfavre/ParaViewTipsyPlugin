set(classes
  vtkTipsyReader)

set(private_headers
  vtkTipsyReader.h)

vtk_module_add_module(TipsyReader
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  Tipsy_Server.xml)

ADD_EXECUTABLE(TestSimpleTipsyReader TestSimpleTipsyReader.cxx)
TARGET_LINK_LIBRARIES(TestSimpleTipsyReader TipsyReader)

