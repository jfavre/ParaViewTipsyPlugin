add_executable(TestSimpleTipsyReader TestSimpleTipsyReader.cxx)
target_link_libraries(TestSimpleTipsyReader
	PUBLIC IOTipsy
	PRIVATE
	  VTK::InteractionStyle
	  VTK::RenderingCore
	  VTK::RenderingOpenGL2
	  )

add_executable(TestTipsyVertexCells TestTipsyVertexCells.cxx)
target_link_libraries(TestTipsyVertexCells
	PUBLIC IOTipsy
	PRIVATE
	  VTK::CommonExecutionModel
	  )

# Tipsy binary file used by TestTipsyVertexCells (local data, machine specific)
set(TIPSY_TEST_FILE "/local/data/Tipsy/hr8799_bol_bd1.017300"
  CACHE FILEPATH "Tipsy binary file used by TestTipsyVertexCells")
add_test(NAME TestTipsyVertexCells
  COMMAND TestTipsyVertexCells -f "${TIPSY_TEST_FILE}")
# the test binary returns 1 (skip) when the data file is not present
set_tests_properties(TestTipsyVertexCells PROPERTIES SKIP_RETURN_CODE 1)

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

