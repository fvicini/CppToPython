cmake_minimum_required(VERSION 3.18)

project(GeDiM4Py LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library test variables
###############################################################################################
set(GEDIMFORPY_SOURCES ${PROJECT_NAME}_sources)
set(GEDIMFORPY_HEADERS ${PROJECT_NAME}_headers)
set(GEDIMFORPY_INCLUDE ${PROJECT_NAME}_include)
set(GEDIMFORPY_LINKED_LIBRARIES ${PROJECT_NAME}_linked_libraries)

# ADD GEDIM DIRECTORY
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gedim/GeDiM)

# Add dependencies
################################################################################

## GeDiM
if (NOT TARGET GeDiM::GeDiM)
	find_package(GeDiM CONFIG REQUIRED)
	list(APPEND THIRDPARTY_LINKED_LIBRARIES GeDiM::GeDiM)
else()
	get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
	get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
	get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
	get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
	get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

	list(APPEND THIRDPARTY_LINKED_LIBRARIES ${GEDIM_LIBRARY_NAME})# ${GEDIM_LIBRARY_LINKED_LIBRARIES})
	list(APPEND THIRDPARTY_HEADERS ${GEDIM_LIBRARY_HEADERS})
	list(APPEND THIRDPARTY_INCLUDE ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
endif()

## Python
find_package(Python3 REQUIRED Development)
list(APPEND THIRDPARTY_LINKED_LIBRARIES Python3::Python)

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)

## Blas
find_package(BLAS REQUIRED)

## Lapack
find_package(LAPACK REQUIRED)

## VTK
if(${ENABLE_VTK})
  find_package(VTK CONFIG REQUIRED)
endif()

## triangle
find_package(triangle CONFIG REQUIRED)

## tetgen
find_package(tetgen CONFIG REQUIRED)

## GTest
if (CMAKE_VERSION VERSION_LESS 3.20)
    find_package(GTest REQUIRED)
    find_package(GMock REQUIRED)

    list(APPEND TEST_LINKED_LIBRARIES GTest::GTest)
    list(APPEND TEST_LINKED_LIBRARIES GTest::Main)
    list(APPEND TEST_LINKED_LIBRARIES ${GMOCK_BOTH_LIBRARIES})
else()
    find_package(GTest REQUIRED)

    list(APPEND TEST_LINKED_LIBRARIES GTest::gtest)
    list(APPEND TEST_LINKED_LIBRARIES GTest::gtest_main)
    list(APPEND TEST_LINKED_LIBRARIES GTest::gmock)
    list(APPEND TEST_LINKED_LIBRARIES GTest::gmock_main)
endif()

# Insert Sources
###############################################################################################

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND ${GEDIMFORPY_HEADERS} ${src_headers})
list(APPEND ${GEDIMFORPY_SOURCES} ${src_sources})
list(APPEND ${GEDIMFORPY_INCLUDE} ${src_include})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)

# Create executable
################################################################################
add_executable(hello_program ${THIRDPARTY_HEADERS} hello.cpp)
target_link_libraries(hello_program ${THIRDPARTY_LINKED_LIBRARIES})
target_include_directories(hello_program SYSTEM PRIVATE ${THIRDPARTY_INCLUDE})

add_library(hello SHARED ${THIRDPARTY_HEADERS} hello.cpp)
target_link_libraries(hello ${THIRDPARTY_LINKED_LIBRARIES})
target_include_directories(hello SYSTEM PRIVATE ${THIRDPARTY_INCLUDE})
target_compile_options(hello PUBLIC -fPIC -Wbad-function-cast -Wcast-function-type)
set_target_properties(
    hello
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "hello"
        LINKER_LANGUAGE CXX
    )

add_executable(GeDiM4Py_Test ${THIRDPARTY_HEADERS} ${${GEDIMFORPY_HEADERS}} ${${GEDIMFORPY_SOURCES}} ${test_sources} GeDiM4Py.cpp)
target_link_libraries(GeDiM4Py_Test ${THIRDPARTY_LINKED_LIBRARIES})
target_link_libraries(GeDiM4Py_Test ${TEST_LINKED_LIBRARIES})
target_include_directories(GeDiM4Py_Test SYSTEM PRIVATE ${THIRDPARTY_INCLUDE})
target_include_directories(GeDiM4Py_Test SYSTEM PRIVATE ${${GEDIMFORPY_INCLUDE}})
target_include_directories(GeDiM4Py_Test SYSTEM PRIVATE ${test_include})

add_library(GeDiM4Py SHARED ${THIRDPARTY_HEADERS} ${${GEDIMFORPY_HEADERS}} ${${GEDIMFORPY_SOURCES}})
target_link_libraries(GeDiM4Py ${THIRDPARTY_LINKED_LIBRARIES})
target_include_directories(GeDiM4Py SYSTEM PRIVATE ${THIRDPARTY_INCLUDE})
target_include_directories(GeDiM4Py SYSTEM PRIVATE ${${GEDIMFORPY_INCLUDE}})
target_compile_options(GeDiM4Py PUBLIC -fPIC)
set_target_properties(
    GeDiM4Py
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "GeDiM4Py"
        LINKER_LANGUAGE CXX
    )
