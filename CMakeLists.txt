cmake_minimum_required(VERSION 2.8)

#===============================================================================
# Regression tests
#===============================================================================

# This allows for dashboard configuration
include(CTest)

# Get a list of all the tests to run
file(GLOB_RECURSE TESTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.py)

# Loop through all the tests
foreach(test ${TESTS})

  # Get test information
  get_filename_component(TEST_NAME ${test} NAME)
  get_filename_component(TEST_PATH ${test} PATH)

  if(${COVERAGE})
    add_test(NAME ${TEST_NAME}
      WORKING_DIRECTORY ${TEST_PATH}
      COMMAND ${PYTHON_EXECUTABLE} run --omit=${OMIT} ${TEST_NAME}
  )
  else()
    add_test(NAME ${TEST_NAME}
      WORKING_DIRECTORY ${TEST_PATH}
      COMMAND ${PYTHON_EXECUTABLE} ${TEST_NAME}
  )
  endif()

endforeach(test)
