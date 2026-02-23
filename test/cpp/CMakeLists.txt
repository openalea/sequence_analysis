if(NOT WITH_TEST)
  return()
endif()

# List all test files
file(GLOB TEST_FILES "*.cpp")

# Loop over test files and create one executable per file
foreach(test_src ${TEST_FILES})
  # Remove the .cpp extension for the target name
  get_filename_component(TEST_NAME ${test_src} NAME_WE)


  add_executable(${TEST_NAME} ${test_src})

  target_include_directories(${TEST_NAME}
    PUBLIC
      ${CONDA_ENV}include
  )

  find_library(STAT_TOOL_LIB NAMES oastat_tool)
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)

  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif("${isSystemDir}" STREQUAL "-1")

  link_directories("${CONDA_ENV}/lib") 
  target_link_directories(${TEST_NAME} PRIVATE "${CONDA_ENV}/lib")
  # Link with your main library and GoogleTest (if used)
  target_link_libraries(${TEST_NAME} PRIVATE oastat_tool)
  target_link_libraries(${TEST_NAME} PRIVATE oasequence_analysis)

  # Register the test with CTest
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

