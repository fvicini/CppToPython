cmake_minimum_required(VERSION 3.18)

project(hello LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library test variables
###############################################################################################
set(HELLO_SOURCES ${PROJECT_NAME}_sources)
set(HELLO_HEADERS ${PROJECT_NAME}_headers)
set(HELLO_INCLUDE ${PROJECT_NAME}_include)
set(HELLO_LINKED_LIBRARIES ${PROJECT_NAME}_linked_libraries)

# ADD GEDIM DIRECTORY
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gedim/GeDiM)

# Add dependencies
################################################################################

## GeDiM
if (NOT TARGET GeDiM::GeDiM)
	find_package(GeDiM CONFIG REQUIRED)
	list(APPEND ${HELLO_LINKED_LIBRARIES} GeDiM::GeDiM)
else()
	get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
	get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
	get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
	get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
	get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

	list(APPEND ${HELLO_LINKED_LIBRARIES} ${GEDIM_LIBRARY_NAME})# ${GEDIM_LIBRARY_LINKED_LIBRARIES})
	list(APPEND ${HELLO_HEADERS} ${GEDIM_LIBRARY_HEADERS})
	list(APPEND ${HELLO_INCLUDE} ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
endif()

## Python
find_package(Python3 REQUIRED Development)
list(APPEND ${HELLO_LINKED_LIBRARIES} Python3::Python)

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

# Insert Sources
###############################################################################################

# ADD SOURCES
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
#list(APPEND ${HELLO_SOURCES} ${src_sources})
list(APPEND ${HELLO_SOURCES} ${${HELLO_HEADERS}})
#list(APPEND ${HELLO_INCLUDE} ${src_include})
list(APPEND ${HELLO_SOURCES} hello.cpp)

# Create executable
################################################################################
add_executable(hello_program ${${HELLO_SOURCES}})
target_link_libraries(hello_program ${${HELLO_LINKED_LIBRARIES}})
target_include_directories(hello_program SYSTEM PRIVATE ${${HELLO_INCLUDE}})

add_library(hello SHARED ${${HELLO_SOURCES}})
target_link_libraries(hello ${${HELLO_LINKED_LIBRARIES}})
target_include_directories(hello SYSTEM PRIVATE ${${HELLO_INCLUDE}})
target_compile_options(hello PUBLIC -fPIC)
set_target_properties(
    hello
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "hello"
        LINKER_LANGUAGE CXX
    )
