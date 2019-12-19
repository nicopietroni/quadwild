cmake_minimum_required(VERSION 2.8.12)
project(cork)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_library(cork STATIC 
  src/file_formats/files.cpp
  src/file_formats/ifs.cpp
  src/file_formats/off.cpp
  src/isct/empty3d.cpp
  src/isct/quantization.cpp
  src/isct/triangle.c
  src/util/log.cpp
  src/util/timer.cpp
  src/cork.cpp)
add_executable(cork-bin src/main.cpp)
set(INCLUDEDIRS
  src/mesh
  src/accel/
  src/file_formats/
  src/isct/
  src/math/
  src/mesh/
  src/rawmesh/
  src/util/
)
target_include_directories(cork PRIVATE ${INCLUDEDIRS})
target_include_directories(cork-bin PRIVATE ${INCLUDEDIRS})
SET_TARGET_PROPERTIES(cork-bin PROPERTIES OUTPUT_NAME cork)
find_package(GMP REQUIRED)
target_link_libraries(cork-bin PRIVATE cork gmp)
target_compile_definitions(cork PRIVATE 
  -DNO_TIMER -DREDUCED -DCDT_ONLY -DTRILIBRARY -DANSI_DECLARATORS )
