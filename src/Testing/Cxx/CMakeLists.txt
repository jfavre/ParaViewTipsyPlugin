add_executable(TestSimpleTipsyReader TestSimpleTipsyReader.cxx)
target_link_libraries(TestSimpleTipsyReader
	PUBLIC IOTipsy
	PRIVATE
	  VTK::InteractionStyle
	  VTK::RenderingCore
	  VTK::RenderingOpenGL2
	  )

