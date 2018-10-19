cmake_minimum_required(VERSION 3.5)

project(TriventProcessor VERSION 0.9.0 LANGUAGES CXX)

# -- Find Dependencies
find_package(ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES)
include(ilcsoft_default_settings)
find_package(Marlin 1.10 REQUIRED) # minimum required Marlin version
find_package(CLHEP REQUIRED) # Required for Marlin >v1.10
find_package(ROOT 6.08 REQUIRED)

# If no nlohmann json installation found, use the embeded files
find_package(nlohmann_json 3.2.0 QUIET )
if (NOT nlohmann_json_FOUND)
  message("-- WARNING: nlohmann_json version>3.2.0 not found, using embeded version")
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  add_subdirectory(nlohmann_json)
endif()

# -- Set some CMake specific variables
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -- Create target and set properties
add_library(TriventProc SHARED EventBuilder/src/TriventProc.cc)
target_include_directories(TriventProc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/EventBuilder/include)

target_compile_features(TriventProc PUBLIC cxx_std_14)

# Define cxxFlags to use
target_compile_options(TriventProc PRIVATE $<$<CXX_COMPILER_ID:Clang>:-Weverything -Wno-c++98-compat -Wno-global-constructors -Wno-exit-time-destructors -Wno-padded -Wdocumentation> -Wall -Wextra -pedantic -Wnon-virtual-dtor -Weffc++ -Wshadow -Wuninitialized -Winit-self -D_FORTIFY_SOURCE)

# ILCSOFT/ROOT uses old style cmake, need to manually set include dirs. (As system to avoid getting crushed by warnings at compilation)
foreach(pkg ROOT Marlin CLHEP)
  list(APPEND TARGET_LIBS ${${pkg}_LIBRARIES})
  target_include_directories(TriventProc SYSTEM PUBLIC ${${pkg}_INCLUDE_DIRS})
endforeach()

target_link_libraries(TriventProc PUBLIC ${TARGET_LIBS} nlohmann_json::nlohmann_json)

# -- Install instruction
install(TARGETS TriventProc
    EXPORT TriventProc-targets
    LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
    ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib
    INCLUDES DESTINATION ${PROJECT_SOURCE_DIR}/include
)

install(EXPORT TriventProc-targets
  FILE
    TriventProcTargets.cmake
  NAMESPACE
    Trivent::
  DESTINATION
    ${PROJECT_SOURCE_DIR}/lib/cmake/TriventProc
)
