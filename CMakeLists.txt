cmake_minimum_required(VERSION 3.18)

project(hello LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add dependencies
################################################################################
find_package(Python3 REQUIRED Development)

# Create executable
################################################################################
add_executable(hello_program hello.cpp)
target_link_libraries(hello_program PUBLIC Python3::Python)

add_library(hello SHARED hello.cpp)
target_link_libraries(hello PUBLIC Python3::Python)
target_compile_options(hello PUBLIC -fPIC)
set_target_properties(
    hello
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "hello"
        LINKER_LANGUAGE CXX
    )
