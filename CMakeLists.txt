cmake_minimum_required(VERSION 3.5)

project(ChaseICING)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

# https://stackoverflow.com/a/6646518/2289030
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost 1.62.0 REQUIRED) # version avail on WSL debian

option(DEBUG_GPROF OFF)
if(${DEBUG_GPROF})
    add_definitions(-pg -no-pie)
    set(CMAKE_EXE_LINKER_FLAGS -pg ${CMAKE_EXE_LINKER_FLAGS})
endif()

include(ExternalProject)

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})

    add_subdirectory(lib)

	add_executable(test_icing test/test_icing.cpp)

    # get the controller
    ExternalProject_Add(imtc
        GIT_REPOSITORY git@gitlab.com:dejournett/imt-c-controller.git
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
    )

    ExternalProject_Get_Property(imtc install_dir)
    include_directories(${install_dir}/src/imtc/include)

    add_dependencies(test_icing imtc)
    target_link_libraries(test_icing ${install_dir}/src/imtc-build/lib/libimt.so)
endif()
