file(GLOB SOURCES "*.cpp")


target_sources(
  _sequence_analysis
  PRIVATE
    ${SOURCES}
)
