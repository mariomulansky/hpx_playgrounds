cmake_minimum_required (VERSION 2.8.4 FATAL_ERROR)
project (test CXX)

set(CMAKE_MODULE_PATH 
       /home/mario/hpx/share/cmake-2.8/Modules
)

# Instruct cmake to find the HPX settings
find_package(HPX)

include_directories( ${HPX_INCLUDE_DIR} )
link_directories( ${HPX_LIBRARY_DIR} )

add_hpx_executable(odeint_dataflow_test
    ESSENTIAL
    SOURCES odeint_dataflow_test.cpp 
    DEPENDENCIES dataflow )

add_hpx_executable(subspace_test
    ESSENTIAL
    SOURCES subspace_test.cpp dataflow_operations.cpp
    DEPENDENCIES dataflow )

add_hpx_executable(subspace_ptr_test
    ESSENTIAL
    SOURCES subspace_ptr_test.cpp 
    DEPENDENCIES dataflow )

add_hpx_executable(dataflow_refs
    ESSENTIAL
    SOURCES dataflow_refs.cpp 
    DEPENDENCIES dataflow )

add_hpx_executable(component_test
    ESSENTIAL
    SOURCES component_test.cpp op.cpp
    )

add_hpx_executable(dataflow_resize
    ESSENTIAL
    SOURCES dataflow_resize.cpp 
    DEPENDENCIES dataflow )
