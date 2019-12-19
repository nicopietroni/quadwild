# This is a minimalist cmake file designed to be used
# with libigl. If you need CoMISo for other purposes
# please use the CMakeLists_original.cmake which
# supports all the features of CoMISo

cmake_minimum_required(VERSION 2.8.12)

project(CoMISo)

add_library(CoMISo STATIC
  Solver/EigenLDLTSolver.hh Solver/EigenLDLTSolver.cc
  Solver/MISolver.hh Solver/MISolver.cc
  Solver/ConstrainedSolver.hh Solver/ConstrainedSolver.cc
  Solver/IterativeSolverT.hh Solver/IterativeSolverT.cc
  Solver/GMM_Tools.hh Solver/GMM_Tools.cc)

target_include_directories(CoMISo PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  "${CMAKE_CURRENT_SOURCE_DIR}/../nanogui/ext/eigen"
  "${CMAKE_CURRENT_SOURCE_DIR}/ext/gmm-4.2/include")
target_compile_definitions(CoMISo PRIVATE -DINCLUDE_TEMPLATES)
if (WIN32)
  target_compile_definitions(CoMISo PRIVATE -D_USE_MATH_DEFINES -DNOMINMAX)
endif(WIN32)
