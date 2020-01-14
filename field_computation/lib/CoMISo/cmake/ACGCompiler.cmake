################################################################################
# Custom settings for compiler flags and similar
################################################################################

if (UNIX)
        
  set (CMAKE_CFLAGS_RELEASE "-O3 -DINCLUDE_TEMPLATES -W -Wall -Wno-unused -DNDEBUG")
  set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DINCLUDE_TEMPLATES -ftemplate-depth-100 -W -Wall -Wno-unused -DNDEBUG")
  set (CMAKE_C_FLAGS_DEBUG "-g -DINCLUDE_TEMPLATES -W -Wall -Wno-unused -DDEBUG")
  set (CMAKE_CXX_FLAGS_DEBUG "-g -DINCLUDE_TEMPLATES -ftemplate-depth-100 -W -Wall -Wno-unused -DDEBUG")

  if (APPLE)
    add_definitions( -DARCH_DARWIN )
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-non-virtual-dtor")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-non-virtual-dtor")
  endif ()
endif ()
