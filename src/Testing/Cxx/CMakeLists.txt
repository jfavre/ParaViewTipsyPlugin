add_executable(TestSimpleTipsyReader TestSimpleTipsyReader.cxx)
target_link_libraries(TestSimpleTipsyReader
	PUBLIC IOTipsy
	PRIVATE
	  VTK::InteractionStyle
	  VTK::RenderingCore
	  VTK::RenderingOpenGL2
	  )

find_package(Viskores REQUIRED QUIET)

add_executable(TestTipsyVTKmConvert TestTipsyVTKmConvert.cxx)
target_link_libraries(TestTipsyVTKmConvert
	PUBLIC IOTipsy
	  viskores::io
	  viskores::filter
	  )

add_executable(TipsyReadAndSwap TipsyReadAndSwap.cxx)
target_link_libraries(TipsyReadAndSwap
	PUBLIC IOTipsy
	  )

