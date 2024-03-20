set(classes
  vtkTipsyReader)

vtk_module_add_module(IOTipsy
  CLASSES ${classes}
  )

paraview_add_server_manager_xmls(
  XMLS  Tipsy_Server.xml)
