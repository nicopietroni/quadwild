cmake_minimum_required (VERSION 2.6)

project(CoMISo)

SET(CMAKE_SKIP_RPATH  TRUE)

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

include (ACGCommon)

if (false)
acg_qt4 ()
# change to 0 if QT should not be used
set( WANT_COMISO_QT 0 )
if( QT4_FOUND)
  #message( WARNING " QT4 FOUND" )
  if( WANT_COMISO_QT )
    add_definitions (-DQT4_FOUND)
  #  message( WARNING " USING QT4" )
  endif ()
  set (COMISO_QT4_CONFIG_FILE_SETTINGS "#define COMISO_QT4_AVAILABLE 1" )
else()
  set (COMISO_QT4_CONFIG_FILE_SETTINGS "#define COMISO_QT4_AVAILABLE 0" )
endif ()
endif(false)

acg_get_version ()

include (ACGOutput)

set(COMISO_INCLUDE_DIRECTORIES "")
set(COMISO_LINK_DIRECTORIES "")
set(COMISO_LINK_LIBRARIES "")
list( APPEND COMISO_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/gmm/include" )
list( APPEND COMISO_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../eigen" )

set (COMISO_Eigen3_CONFIG_FILE_SETTINGS "#define COMISO_Eigen3_AVAILABLE 1" )
set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 1" )


# We require cgal with its blas on windows
if (WIN32)
  find_package(CGAL)
  if (CGAL_FOUND)
    list( APPEND COMISO_INCLUDE_DIRECTORIES ${CGAL_INCLUDE_DIR} )
    list( APPEND COMISO_LINK_DIRECTORIES   ${CGAL_LIBRARY_DIR} )
    list( APPEND COMISO_LINK_LIBRARIES  ${CGAL_LIBRARIES} )
  endif()
endif()

find_package (BLAS)
if (BLAS_FOUND )
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 1" )
  
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_DIRECTORIES   ${BLAS_LIBRARY_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${BLAS_LIBRARIES} )
else()
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 0" )
  message (FATAL_ERROR "BLAS not found!")
endif ()

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  ${CMAKE_CURRENT_BINARY_DIR}
  ${COMISO_INCLUDE_DIRECTORIES}
)

# generate dllexport macros on windows
if (WIN32)
  add_definitions(-DCOMISODLL)
endif ()


link_directories (
  ${COMISO_LINK_DIRECTORIES}
)

# source code directories
set (directories 
  .
  Solver
  NSolver
  EigenSolver
  Config
  Utils
  QtWidgets
)

# generate dllexport macros on windows
if (WIN32)
  add_definitions(-DCOMISODLL)
endif ()

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )
endif ()

# collect all header,source and ui files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})
acg_append_files (ui "*.ui" ${directories})


macro (of_list_filter _list)
  if (WIN32)
   foreach (_element ${${_list}})
     if (_element MATCHES "gnuplot_i\\.(cc|hh)$")
       list (REMOVE_ITEM ${_list} ${_element})
     endif ()
   endforeach ()
  endif ()
endmacro ()

of_list_filter  ( headers )
of_list_filter  ( sources )


# remove template cc files from source file list
acg_drop_templates (sources)

if( QT4_FOUND)
# genereate uic and moc targets
acg_qt4_autouic (uic_targets ${ui})
acg_qt4_automoc (moc_targets ${headers})
endif()
acg_add_library (CoMISo SHARED ${uic_targets} ${sources} ${headers} ${moc_targets})

if (NOT APPLE)
	target_link_libraries (CoMISo 
	  ${QT_LIBRARIES}
	  ${COMISO_LINK_LIBRARIES}
	)
else(NOT APPLE)
	target_link_libraries (CoMISo
	  ${QT_LIBRARIES}
	  ${COMISO_LINK_LIBRARIES}
	)
endif(NOT APPLE)

# display results
acg_print_configure_header (COMISO "CoMISo")

# write config file
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/Config/config.hh.in"
     "${CMAKE_CURRENT_SOURCE_DIR}/Config/config.hh" @ONLY IMMEDIATE)
     


