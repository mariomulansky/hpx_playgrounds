cmake_minimum_required (VERSION 2.8.4 FATAL_ERROR)
project (test CXX)

set(CMAKE_MODULE_PATH 
                      $ENV{HPX_ROOT}/share/cmake-2.8/Modules
)

# Instruct cmake to find the HPX settings
find_package(HPX)

include_directories( /home/mario/odeint-v2 )
include_directories( ${HPX_INCLUDE_DIR} )
link_directories( ${HPX_LIBRARY_DIR} )

add_hpx_executable(test_local_dataflow
    ESSENTIAL
    SOURCES test_local_dataflow.cpp
    DEPENDENCIES iostreams
)

add_hpx_executable(test2d
    ESSENTIAL
    SOURCES test2d.cpp
    DEPENDENCIES iostreams
)

add_hpx_executable(perf2d
    ESSENTIAL
    SOURCES perf2d.cpp
    DEPENDENCIES iostreams
)

add_hpx_executable(const_test
    ESSENTIAL
    SOURCES const_test.cpp
    DEPENDENCIES iostreams
)

add_hpx_executable(test_const_functor
    ESSENTIAL
    SOURCES test_const_functor.cpp
    DEPENDENCIES iostreams
)

add_hpx_executable(test_shared_dataflow
    ESSENTIAL
    SOURCES test_shared_dataflow.cpp
    DEPENDENCIES iostreams
)

add_hpx_executable(test_swap
    ESSENTIAL
    SOURCES test_swap.cpp
    DEPENDENCIES iostreams
)
