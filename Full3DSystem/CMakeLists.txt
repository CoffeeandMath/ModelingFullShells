##
#  CMake script for the step-8 tutorial program:
##
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE '-O3')

add_definitions(-w)


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
project(CylindricalShells)
FIND_PACKAGE(deal.II 9.2.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/include
)



ADD_LIBRARY(elasticproblemclass ${CMAKE_SOURCE_DIR}/src/elastic_problem.cpp)
DEAL_II_SETUP_TARGET(elasticproblemclass)

ADD_LIBRARY(read_file ${CMAKE_SOURCE_DIR}/src/read_file.cpp)
DEAL_II_SETUP_TARGET(read_file)


ADD_EXECUTABLE(run ${CMAKE_SOURCE_DIR}/src/main.cpp)
DEAL_II_SETUP_TARGET(run)
TARGET_LINK_LIBRARIES(elasticproblemclass read_file)
TARGET_LINK_LIBRARIES(run elasticproblemclass)

