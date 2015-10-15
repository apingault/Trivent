# Install script for directory: /Users/antoine/CALICE/DataAnalysis/triventArnaud/EventBuilder

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/Users/antoine/CALICE/DataAnalysis/triventArnaud")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE DIRECTORY FILES "/Users/antoine/CALICE/DataAnalysis/triventArnaud/EventBuilder/./include" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\~$" EXCLUDE REGEX "/[^/]*\\#[^/]*$" EXCLUDE REGEX "/\\.\\#[^/]*$" EXCLUDE REGEX "/[^/]*cvs$" EXCLUDE REGEX "/[^/]*\\.svn$" EXCLUDE)
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES
    "/Users/antoine/CALICE/DataAnalysis/triventArnaud/lib/libTrivent.0.4.1.dylib"
    "/Users/antoine/CALICE/DataAnalysis/triventArnaud/lib/libTrivent.0.4.dylib"
    "/Users/antoine/CALICE/DataAnalysis/triventArnaud/lib/libTrivent.dylib"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTrivent.0.4.1.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTrivent.0.4.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTrivent.dylib"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
        -id "/Users/antoine/CALICE/DataAnalysis/triventArnaud/lib/libTrivent.0.4.dylib"
        "${file}")
      execute_process(COMMAND /usr/bin/install_name_tool
        -add_rpath "/Users/antoine/CALICE/DataAnalysis/triventArnaud/lib"
        "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

