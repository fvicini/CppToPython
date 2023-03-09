list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/GeDiM4Py_Logic.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/GeDiM4Py_Interface.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/FEM_RefElement_Langrange_PCC_Triangle_2D.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/PDE_Equation.hpp)

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/GeDiM4Py_Logic.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/GeDiM4Py_Interface.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/FEM_RefElement_Langrange_PCC_Triangle_2D.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/PDE_Equation.cpp)

list(APPEND src_include ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_include ${src_include} PARENT_SCOPE)
