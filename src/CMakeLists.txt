# @author  Hector Rieiro


set(APP_HEADERS_PROPS
"${CMAKE_CURRENT_BINARY_DIR}/../include/Analysis/Version.h"
../include/Analysis/Numerics.h
../include/Analysis/FixationAnalysis.h
../include/Analysis/MicrosaccadeAnalysis.h
../include/Analysis/PursuitAnalysis.h
../include/Analysis/FixationEstabilityAnalysis.h
../include/Analysis/GuidedSaccadeAnalysis.h
../include/Analysis/PupillometryAnalysis.h
../include/Analysis/SmoothConvergenceAnalysis.h
../include/Analysis/OKNAnalysis.h
../include/Analysis/AnalysisImpl.h
../include/Analysis/AnalysisExport.h
../include/Analysis/DataTypes.h
../include/Analysis/DataFile.h
../include/Analysis/HDF5.h
)

set(APP_HEADERS
 ../include/Analysis/Version.h.in
../include/Analysis/Numerics.h
../include/Analysis/Algorithms.h
../include/Analysis/FixationAnalysis.h
../include/Analysis/MicrosaccadeAnalysis.h
../include/Analysis/PursuitAnalysis.h
../include/Analysis/GuidedSaccadeAnalysis.h
../include/Analysis/FixationEstabilityAnalysis.h
../include/Analysis/GuidedSaccadeAnalysis.h
../include/Analysis/PupillometryAnalysis.h
../include/Analysis/OKNAnalysis.h
../include/Analysis/SmoothConvergenceAnalysis.h
../include/Analysis/AnalysisImpl.h
../include/Analysis/AnalysisExport.h
../include/Analysis/DataTypes.h
../include/Analysis/DataFile.h
../include/Analysis/HDF5.h
)

set(APP_SOURCES
  Numerics.cpp Algorithms.cpp FixationAnalysis.cpp MicrosaccadeAnalysis.cpp  
  FixationEstabilityAnalysis.cpp PursuitAnalysis.cpp GuidedSaccadeAnalysis.cpp PupillometryAnalysis.cpp SmoothConvergenceAnalysis.cpp OKNAnalysis.cpp DataFile.cpp HDF5.cpp
  AnalysisExport.cpp)

find_package(HDF5 NAMES hdf5 COMPONENTS CXX static)

set(LIBRARY_NAME Analysis)
add_library(${LIBRARY_NAME} SHARED ${APP_SOURCES} ${APP_HEADERS})
target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include/" "${CMAKE_CURRENT_BINARY_DIR}/../include/" ${Boost_INCLUDE_DIRS}  ${HDF5_INCLUDE_DIR} )
# Links the target exe (projBrunner) against the library (libprojB) and other required libraries
target_link_libraries(${LIBRARY_NAME} PUBLIC ${Boost_LIBRARIES} jsoncpp_lib_static ceres hdf5_cpp-static zlibstatic)
set_property(TARGET ${LIBRARY_NAME} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

add_executable(TestHDF5Conversion TestHDF5Conversion.cpp)
target_link_libraries(TestHDF5Conversion Analysis  ${Boost_LIBRARIES})
target_include_directories(TestHDF5Conversion PUBLIC ${Boost_INCLUDE_DIRS})
set_property(TARGET TestHDF5Conversion PROPERTY FOLDER "Tests")

add_executable(TestSegmentation TestSegmentation.cpp)
target_link_libraries(TestSegmentation Analysis  ${Boost_LIBRARIES})
target_include_directories(TestSegmentation PUBLIC ${Boost_INCLUDE_DIRS})
set_property(TARGET TestSegmentation PROPERTY FOLDER "Tests")

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qpar /Qpar-report:1 /Qvec-report:1")
endif(MSVC)

# Says how and where to install software
# Targets:
#   * <prefix>/lib/<libraries>
#   * header location after install: <prefix>/include/<project>/*.h
#   * headers can be included by C++ code `#<project>/Bar.hpp>`
#install(TARGETS ${LIBRARY_NAME} ${APPLICATION_NAME}
#  EXPORT ${PROJECT_NAME}Targets            # for downstream dependencies
#  ARCHIVE DESTINATION lib COMPONENT libs   # static lib
#  RUNTIME DESTINATION bin COMPONENT libs   # binaries
#  LIBRARY DESTINATION lib COMPONENT libs   # shared lib
#  FRAMEWORK DESTINATION bin COMPONENT libs # for mac
#  PUBLIC_HEADER DESTINATION include/${PROJECT_NAME} COMPONENT devel   # headers for mac (note the different component -> different package)
#  INCLUDES DESTINATION include             # headers
#)
