# - Find Eigen3
# Find the native GMM headers and libraries.
#
#  Eigen3_INCLUDE_DIR -  where to find <Eigen/Dense>, etc.
#  Eigen3_FOUND        - True if Eigen3 found.

IF (Eigen3_INCLUDE_DIR)
  # Already in cache, be silent
  SET(Eigen3_FIND_QUIETLY TRUE)
ENDIF (Eigen3_INCLUDE_DIR)

GET_FILENAME_COMPONENT(module_file_path ${CMAKE_CURRENT_LIST_FILE} PATH )

# Look for the header file.
FIND_PATH(Eigen3_INCLUDE_DIR NAMES Eigen/Dense 
                             PATHS /usr/include/eigen3
                                   /usr/local/include
                                   /usr/local/include/eigen3/
                                   /opt/local/include/eigen3/
                                   "c:\\libs\\eigen3\\include"
		                   "c:\\libs\\eigen\\include"
				    ${PROJECT_SOURCE_DIR}/MacOS/Libs/eigen3/include
                                   ../../External/include
                                   ${module_file_path}/../../../External/include)


# Copy the results to the output variables.
IF(Eigen3_INCLUDE_DIR )
  SET(Eigen3_FOUND 1)
  SET(Eigen3_INCLUDE_DIR ${Eigen3_INCLUDE_DIR})
ELSE(Eigen3_INCLUDE_DIR )
  SET(Eigen3_FOUND 0)
  SET(Eigen3_INCLUDE_DIR)
ENDIF(Eigen3_INCLUDE_DIR )

# Report the results.
IF(NOT Eigen3_FOUND)
  SET(Eigen3_DIR_MESSAGE
    "Eigen3 was not found. Make sure Eigen3_INCLUDE_DIR is set to the directories containing the include and lib files for Eigen3. .")
ELSE (NOT Eigen3_FOUND)
  IF(NOT Eigen3_FIND_QUIETLY)
    MESSAGE(STATUS "Looking for Eigen3 - found")
  ENDIF(NOT Eigen3_FIND_QUIETLY)
ENDIF(NOT Eigen3_FOUND)

