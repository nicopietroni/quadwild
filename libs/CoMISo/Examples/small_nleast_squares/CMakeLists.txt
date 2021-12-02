include (ACGCommon)
include (CoMISoExample)

# source code directories
set (directories 
  . 
)

# collect all header and source files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})

# remove template cc files from source file list
acg_drop_templates (sources)

if (WIN32)
  acg_add_executable (small_nleast_squares WIN32 ${sources} ${headers} )
elseif (APPLE)
  # generate bundle on mac
  acg_add_executable (small_nleast_squares MACOSX_BUNDLE ${sources} ${headers} )
else ()
  acg_add_executable (small_nleast_squares ${sources} ${headers} )
endif ()

# enable rpath linking
set_target_properties(small_nleast_squares PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

target_link_libraries (small_nleast_squares
  CoMISo
  ${COMISO_LINK_LIBRARIES}
)

if (APPLE)
   # create bundle in "Build" directory and set icon
   # no install needed here, because the whole bundle will be installed in the
   # toplevel CMakeLists.txt
   set_target_properties (
      small_nsolver PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build"
      MACOSX_BUNDLE_INFO_STRING "CoMISo small_nleast_squares"
   )
endif ()
