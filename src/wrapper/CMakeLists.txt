file(GLOB SOURCES "*.cpp")

add_library(_sequence_analysis MODULE ${SOURCES})

target_include_directories(_sequence_analysis PRIVATE
  ${CMAKE_SOURCE_DIR}/src/cpp
  ${SKLBUILD_HEADERS_DIR}/../openalea.stat_tool
  ${Boost_INCLUDE_DIRS}
  ${Python3_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# cf. https://scikit-build-core.readthedocs.io/en/latest/guide/dynamic_link.html#manual-patching
if(APPLE)
  set(origin_token "@loader_path")
else()
  set(origin_token "$ORIGIN")
endif()

# Link to vpsequence_analysis and Boost.Python
target_link_libraries(_sequence_analysis vpsequence_analysis)
target_compile_definitions(_sequence_analysis PRIVATE BOOST_ALL_NO_LIB)
target_link_libraries(_sequence_analysis ${Boost_LIBRARIES})
if(APPLE)
  target_link_libraries(_sequence_analysis "-undefined dynamic_lookup")
else()
  target_link_libraries(_sequence_analysis Python3::Python)
endif()

set_target_properties(_sequence_analysis PROPERTIES
  INSTALL_RPATH "${origin_token}"
)

add_dependencies(_sequence_analysis vpsequence_analysis)
## wrapper_install(_sequence_analysis sequence_analysis)
message(STATUS "Installing wrappers to ${SKBUILD_PLATLIB_DIR}/openalea/sequence_analysis")
set_target_properties(_sequence_analysis PROPERTIES PREFIX "")
install(TARGETS _sequence_analysis
  LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/openalea/sequence_analysis"
)
