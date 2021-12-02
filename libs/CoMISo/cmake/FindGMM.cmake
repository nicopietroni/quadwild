# - Find GMM
# Find the native GMM headers and libraries.
#
#  GMM_INCLUDE_DIR -  where to find GMM.h, etc.
#  GMM_LIBRARIES    - List of libraries when using GMM.
#  GMM_FOUND        - True if GMM found.

IF (GMM_INCLUDE_DIR)
  # Already in cache, be silent
  SET(GMM_FIND_QUIETLY TRUE)
ENDIF (GMM_INCLUDE_DIR)

GET_FILENAME_COMPONENT(module_file_path ${CMAKE_CURRENT_LIST_FILE} PATH )

# Look for the header file.
FIND_PATH(GMM_INCLUDE_DIR NAMES gmm/gmm.h 
                           PATHS /usr/include
                                 /usr/include
                                 /usr/local/include
                                 ~/sw/gmm-4.1/include
				 "c:\\libs\\gmm-4.1\\include"
                                 "c:\\libs\\gmm-3.0\\include"
			         ${PROJECT_SOURCE_DIR}/MacOS/Libs/gmm-3.1/include
                                 ../../External/include
                                 /ACG/acgdev/gcc-4.3-i686/gmm-4.1/include/
                                 /ACG/acgdev/gcc-4.3-i686/gmm-3/include/
                                 ${module_file_path}/../../../External/include)


# Copy the results to the output variables.
IF(GMM_INCLUDE_DIR )
  SET(GMM_FOUND 1)
  SET(GMM_INCLUDE_DIR ${GMM_INCLUDE_DIR})
  IF (WIN32)
  add_definitions(-D_SCL_SECURE_NO_DEPRECATE)
  ENDIF(WIN32)
ELSE(GMM_INCLUDE_DIR )
  SET(GMM_FOUND 0)
  SET(GMM_INCLUDE_DIR)
ENDIF(GMM_INCLUDE_DIR )

# Report the results.
IF(NOT GMM_FOUND)
  SET(GMM_DIR_MESSAGE
    "GMM was not found. Make sure GMM_INCLUDE_DIR is set to the directories containing the include and lib files for GMM. .")
  IF(GMM_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "${GMM_DIR_MESSAGE}")
  ELSEIF(NOT GMM_FIND_QUIETLY)
    MESSAGE(STATUS "${GMM_DIR_MESSAGE}")
  ELSE(NOT GMM_FIND_QUIETLY)
  ENDIF(GMM_FIND_REQUIRED)
ELSE (NOT GMM_FOUND)
  IF(NOT GMM_FIND_QUIETLY)
    MESSAGE(STATUS "Looking for GMM - found")
  ENDIF(NOT GMM_FIND_QUIETLY)
ENDIF(NOT GMM_FOUND)

