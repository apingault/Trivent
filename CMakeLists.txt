cmake_minimum_required(VERSION 3.5)

project(TriventProcessor VERSION 0.9.0 LANGUAGES CXX)
set(TARGET_NAME "TriventProc")

# -- Find Dependencies
find_package(ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES)
include(ilcsoft_default_settings)
find_package(Marlin 1.10 REQUIRED) # minimum required Marlin version
find_package(CLHEP REQUIRED) # Required for Marlin >v1.10
find_package(ROOT 6.08 REQUIRED)

# If no nlohmann json installation found, use the embeded files
find_package(nlohmann_json 3.1.2 QUIET )
if (NOT nlohmann_json_FOUND)
  message("-- WARNING: nlohmann_json version>3.1.2 not found, using embeded version")
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  add_subdirectory(nlohmann_json)
endif()

# -- Set some CMake specific variables
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -- Create target and set properties
add_library(${TARGET_NAME} SHARED EventBuilder/src/TriventProc.cc)
target_include_directories(${TARGET_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/EventBuilder/include>
)

set_target_properties(${TARGET_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

if (${CMAKE_VERSION} VERSION_GREATER 3.8) # cvmfs only has CMAKE 3.6 as of today...
  target_compile_features(${TARGET_NAME} PRIVATE cxx_std_14)
else()
  set_property(TARGET ${TARGET_NAME} PROPERTY CXX_STANDARD 14)
endif()

# Define cxxFlags to use
# Cannot use the sanitize lags with clang as it needs to build the lib as static and marlin needs a shared lib to function...
target_compile_options(${TARGET_NAME} PRIVATE -Wall -Wextra -pedantic -Wnon-virtual-dtor -Weffc++ -Wshadow -Wuninitialized -Winit-self -D_FORTIFY_SOURCE -D_GLIBCXX_ASSERTIONS -fasynchronous-unwind-tables -fexceptions -fpie -fpic
    $<$<CXX_COMPILER_ID:Clang>:-Wno-c++98-compat -Wno-global-constructors -Wno-exit-time-destructors -Wno-padded -Wno-documentation>
    # $<$<CXX_COMPILER_ID:GNU>:-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=undefined>
)

# ILCSOFT/ROOT uses old style cmake, need to manually set include dirs. (As system to avoid getting crushed by warnings at compilation)
foreach(pkg ROOT Marlin CLHEP)
  list(APPEND TARGET_LIBS ${${pkg}_LIBRARIES})
  target_include_directories(${TARGET_NAME} SYSTEM PUBLIC ${${pkg}_INCLUDE_DIRS})
endforeach()

if (DEFINED ${nlohmann_json_VERSION} AND ${nlohmann_json_VERSION} VERSION_GREATER 3.1.2)
  target_link_libraries(${TARGET_NAME} PUBLIC ${TARGET_LIBS} nlohmann_json::nlohmann_json -Wl,-pie)
else()
  target_link_libraries(${TARGET_NAME} PUBLIC ${TARGET_LIBS} nlohmann_json -Wl,-pie)
endif()

# -- Install instruction
install(TARGETS ${TARGET_NAME}
    EXPORT ${TARGET_NAME}-targets
    LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib/marlin_dll
    ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib
    # INCLUDES DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include
)

install(EXPORT ${TARGET_NAME}-targets
  FILE
    ${TARGET_NAME}Targets.cmake
  NAMESPACE
    Trivent::
  DESTINATION
    ${PROJECT_SOURCE_DIR}/lib/cmake/${TARGET_NAME}
)
