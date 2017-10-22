########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
########################################################

# project name
PROJECT(TriventProcessor)

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 5 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# Require C++11
INCLUDE( CheckCXXCompilerFlag )
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
IF( COMPILER_SUPPORTS_CXX11 )
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
  MESSAGE( SEND_ERROR " ${PROJECT_NAME} Requires C++11 support. Please upgrade your compiler !" )
ENDIF()

# Add Rpath for osx
SET( CMAKE_MACOSX_RPATH ON )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.10 REQUIRED ) # minimum required Marlin version
FIND_PACKAGE( CLHEP REQUIRED ) # Required for Marlin >v1.10
FIND_PACKAGE( ROOT 6.08 REQUIRED )

FOREACH( pkg ROOT Marlin CLHEP)
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -O2 -Wno-effc++ -Wno-shadow")

ADD_SUBDIRECTORY( ./EventBuilder )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
