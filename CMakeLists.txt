cmake_minimum_required(VERSION 3.5)

project(ChaseICING)

set(CMAKE_CXX_STANDARD 14)

# https://stackoverflow.com/a/6646518/2289030
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost 1.62.0) # version avail on WSL debian

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	add_executable(test_icing test/test_icing.cpp)
endif()
