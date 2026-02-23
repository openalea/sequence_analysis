# Glob includes and sources
file(GLOB HEADERS "*.h" "*.hpp")
file(GLOB SOURCES "*.cpp")


target_include_directories(oasequence_analysis PUBLIC .)

target_sources(
  oasequence_analysis
  PRIVATE
    ${SOURCES}
)

# Add the headers to be installed with the library

# Add the headers to be installed with the library
install(
  FILES ${HEADERS}
  DESTINATION "${CONDA_ENV}include/sequence_analysis"
)