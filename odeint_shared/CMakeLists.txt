cmake_minimum_required (VERSION 2.8.4 FATAL_ERROR)
project (test CXX)

set(CMAKE_MODULE_PATH 
       /home/mario/hpx/share/cmake-2.8/Modules
)

# Instruct cmake to find the HPX settings
find_package(HPX)

include_directories( ${HPX_INCLUDE_DIR} )
link_directories( ${HPX_LIBRARY_DIR} )

add_hpx_executable(test
    ESSENTIAL
    SOURCES test.cpp dataflow_shared_operations.cpp
    DEPENDENCIES dataflow )

add_hpx_executable(test_shared_dataflow
    ESSENTIAL
    SOURCES test_shared_dataflow.cpp
    DEPENDENCIES dataflow )
