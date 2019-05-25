cmake_minimum_required(VERSION 3.8)

project(ChaseICING)

set(CMAKE_CXX_STANDARD 17)

# https://stackoverflow.com/a/6646518/2289030
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost 1.70.0 COMPONENTS numeric)

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	add_executable(test_icing test/test_icing.cpp)
	target_link_libraries(test_icing ${Boost_LIBRARIES})
endif()
