add_executable(TestSimpleTipsyReader TestSimpleTipsyReader.cxx)
target_link_libraries(TestSimpleTipsyReader
	PUBLIC IOTipsy
	PRIVATE
	  VTK::InteractionStyle
	  VTK::RenderingCore
	  VTK::RenderingOpenGL2
	  )

find_package(VTKm REQUIRED QUIET)

add_executable(TestTipsyVTKmConvert TestTipsyVTKmConvert.cxx)
target_link_libraries(TestTipsyVTKmConvert
	PUBLIC IOTipsy
	  vtkm::io
	  vtkm::filter
	  )

