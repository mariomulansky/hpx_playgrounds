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

add_hpx_executable(test2d
    ESSENTIAL
    SOURCES test2d.cpp
    DEPENDENCIES dataflow iostreams
)

add_hpx_executable(perf2d
    ESSENTIAL
    SOURCES perf2d.cpp
    DEPENDENCIES dataflow iostreams
)
