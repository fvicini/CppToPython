list(APPEND test_sources ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry.hpp)
list(APPEND test_sources ${CMAKE_CURRENT_SOURCE_DIR}/test_python.hpp)

list(APPEND test_include ${CMAKE_CURRENT_SOURCE_DIR})

set(test_sources ${test_sources} PARENT_SCOPE)
set(test_include ${test_include} PARENT_SCOPE)
