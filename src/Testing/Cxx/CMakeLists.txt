ADD_EXECUTABLE(TestSimpleTipsyReader TestSimpleTipsyReader.cxx)
TARGET_LINK_LIBRARIES(TestSimpleTipsyReader
	PUBLIC TipsyReader
	PRIVATE
	  VTK::InteractionStyle
	  VTK::RenderingCore
	  VTK::RenderingOpenGL2
	  )

