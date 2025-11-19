# Glob includes and sources
file(GLOB HEADERS "*.h" "*.hpp")
file(GLOB SOURCES "*.cpp")

add_library(vpsequence_analysis SHARED ${SOURCES})

target_compile_definitions(vpsequence_analysis PRIVATE MESSAGE)

# Optionally add efence
if(WITH_EFENCE)
  target_compile_definitions(vpsequence_analysis PRIVATE DEBUG)
  target_link_libraries(vpsequence_analysis PRIVATE efence)
endif()

message(STATUS "SKBUILD_HEADERS_DIR: ${SKBUILD_HEADERS_DIR}")
target_include_directories(vpsequence_analysis PUBLIC
  ${SKBUILD_HEADERS_DIR}/../openalea.stat_tool
  $ENV{CONDA_PREFIX}/include/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/openalea.stat_tool
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

set(SEQUENCE_ANALYSIS_LIBDIR "${SKBUILD_PLATLIB_DIR}/openalea/sequence_analysis")

message(STATUS "Installing vpsequence_analysis to ${SEQUENCE_ANALYSIS_LIBDIR}")
install(TARGETS vpsequence_analysis
  LIBRARY DESTINATION "${SEQUENCE_ANALYSIS_LIBDIR}"
)

message(STATUS "Installing headers: ${HEADERS}")
message(STATUS "Installing headers to ${SKBUILD_HEADERS_DIR}/sequence_analysis")
install(FILES ${HEADERS}
  DESTINATION "${SKBUILD_HEADERS_DIR}/sequence_analysis"
)
