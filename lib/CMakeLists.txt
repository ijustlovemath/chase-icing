add_library(xyparser STATIC XYParser.cpp)
add_library(interpolator STATIC Interpolator.cpp)
add_library(inversecdf STATIC InverseCDFProcess.cpp)

#target_link_libraries(interpolator xyparser)
#target_link_libraries(inversecdf interpolator xyparser)
