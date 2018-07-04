# @author Hector Rieiro

cmake_minimum_required(VERSION 3.0)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project(SOAT CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake) # main (top) cmake dir
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) # project specific cmake dir
set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?
set(CXX_STANDARD_REQUIRED)

# CMake useful variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib") 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)


find_package(Boost REQUIRED COMPONENTS filesystem system iostreams)
find_package(Threads REQUIRED)


# Set the version number of your project here (format is MAJOR.MINOR.PATCHLEVEL - e.g. 1.0.0)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

include(Common)

# SVN (and its revision)
#include(FindSubversion)
#IF(Subversion_FOUND)
#    if(EXISTS "${CMAKE_SOURCE_DIR}/.svn")
#        Subversion_WC_INFO(${CMAKE_SOURCE_DIR} MY)
#        SET(VCS_REVISION "${MY_WC_REVISION}")
#    else ()
#        SET(VCS_REVISION "-1")
#    endif()
#ELSE(Subversion_FOUND)
#    SET(VCS_REVISION "-1")
#ENDIF(Subversion_FOUND)

## Git (and its revision)
find_package(Git QUIET) # if we don't find git or FindGit.cmake is not on the system we ignore it.
## GetGitRevisionDescription module to retreive branch and revision information from Git
## Starting with Git 1.9 the module will be part of official cMake distribution, until then it has to be
## part of the application
## The Git module will trigger a reconfiguration for each pull that will bring a new revision on the local repository
set (VCS_REVISION "-1")
if(GIT_FOUND)
    include(GetGitRevisionDescription)
    get_git_head_revision(GIT_REFSPEC GIT_SHA1)
    message(STATUS "GIT branch ${GIT_REFSPEC}")
    message(STATUS "GIT revision ${GIT_SHA1}")
    set (VCS_REVISION ${GIT_SHA1})
endif()

# This line finds the boost lib and headers. 
#set(Boost_NO_BOOST_CMAKE ON) #  Don't do a find_package in config mode before searching for a regular boost install.
#find_package(Boost COMPONENTS unit_test_framework program_options REQUIRED)

# This sets the include directory for the reference project. This is the -I flag in gcc.
# All the includes should be in this variable To add new directories to the include path,
# just add them in this list (after a space), or add an include_directories statement in a
# subdirectory's CMakeLists.
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include   # add the binary tree to the search path so that we will find Version.h
)


# Mac needed variables (adapt according to your needs)
set(CMAKE_MACOSX_RPATH ON)
#set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if (NOT TARGET gtest)
	# Enable ExternalProject CMake module
	include(ExternalProject)
	
	# Download and unpack googletest at configure time
	configure_file(googletest-download.in googletest/download/CMakeLists.txt)
	execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
		WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest/download" )
	execute_process(COMMAND "${CMAKE_COMMAND}" --build .
		WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest/download" )
	
	# Prevent GoogleTest from overriding our compiler/linker options
	# when building with Visual Studio
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	
	# Add googletest directly to our build. This adds
	# the following targets: gtest, gtest_main, gmock
	# and gmock_main
	add_subdirectory("${CMAKE_BINARY_DIR}/googletest/src"
					"${CMAKE_BINARY_DIR}/googletest/build")
	
	# The gtest/gmock targets carry header search path
	# dependencies automatically when using CMake 2.8.11 or
	# later. Otherwise we have to add them here ourselves.
	if(CMAKE_VERSION VERSION_LESS 2.8.11)
		include_directories("${gtest_SOURCE_DIR}/include"
							"${gmock_SOURCE_DIR}/include")
	endif()
	set_property(TARGET gtest PROPERTY FOLDER "External")
	set_property(TARGET gtest_main PROPERTY FOLDER "External")
	set_property(TARGET gmock PROPERTY FOLDER "External")
	set_property(TARGET gmock_main PROPERTY FOLDER "External")
endif()

enable_testing()

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(ext/ceres-solver)
set_property(TARGET ceres PROPERTY FOLDER "External")
add_subdirectory(test)
add_subdirectory(ext/jsoncpp)
set_property(TARGET jsoncpp_check PROPERTY FOLDER "External")
set_property(TARGET jsoncpp_lib_static PROPERTY FOLDER "External")
set_property(TARGET jsoncpp_readerwriter_tests PROPERTY FOLDER "External")
set_property(TARGET jsoncpp_test PROPERTY FOLDER "External")
set_property(TARGET jsontestrunner_exe PROPERTY FOLDER "External")
add_subdirectory(ext/zlib-1.2.11)
include_directories("${ZLIB_INCLUDE_DIRS}")
set(ZLIB_LIBRARIES zlibstatic)
set(TARGET zlib PROPERTY FOLDER "External/Zlib")
set(TARGET zlibstatic PROPERTY FOLDER "External/Zlib")

set(HDF5_EXTERNALLY_CONFIGURED 1)
if (HDF5_ENABLE_Z_LIB_SUPPORT)
	set(HDF5_LIB_DEPENDENCIES zlibstatic)
	set(H5_ZLIB_HEADER "zlib.h")
endif()
add_subdirectory(ext/hdf5-1.10.2)
set(TARGET PROPERTY hdf5-static FOLDER "External/HDF5")
set(TARGET PROPERTY hdf5_hl_cpp-static FOLDER "External/HDF5")
set(TARGET PROPERTY hdf5_hl-static FOLDER "External/HDF5")
set(TARGET PROPERTY H5detect FOLDER "External/HDF5")
set(TARGET PROPERTY H5make_libsettings FOLDER "External/HDF5")
set(TARGET PROPERTY hdf5_cpp-static FOLDER "External/HDF5")
#set_global_variable (HDF5_LIBRARIES_TO_EXPORT "hdf5_cpp-static")

#
# SYSTEM PACKAGING (RPM, TGZ, ...)
# _____________________________________________________________________________

include(CPackConfig)

#
# CMAKE PACKAGING (for other CMake projects to use this one easily)
# _____________________________________________________________________________

include(PackageConfigurator)
