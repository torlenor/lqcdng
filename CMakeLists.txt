cmake_minimum_required(VERSION 2.8)

project( LQCDNG ) # this is case sensitive 

######################################

include_directories (${LQCDNG_SOURCE_DIR})
link_directories (${LQCDNG_BINARY_DIR})

SET(CODE_FILES_SU2 puresu2gauge.cc globalsettings.cc su2.cc init.cc simulationkernels/genericsimclass.cc simulationkernels/puresu2gauge.cc measuretime.cc)
SET(CODE_FILES_SU3 puresu3gauge.cc globalsettings.cc su3.cc init.cc simulationkernels/genericsimclass.cc simulationkernels/puresu3gauge.cc measuretime.cc)

set(CMAKE_BUILD_TYPE Release)

SET(GCC_ADDITIONAL_COMPILE_FLAGS  "-O2 -Wall -std=c++11 -DDEBUG")
# SET(GCC_ADDITIONAL_COMPILE_FLAGS  "-O3 -Wall -funroll-loops -std=c++11 -DNDEBUG")
# SET(GCC_ADDITIONAL_COMPILE_FLAGS  "-O3 -Wall -funroll-loops -std=c++11 -DNDEBUG -DPREFETCH")
SET(GCC_ADDITIONAL_LINK_FLAGS     "")
SET(GCC_ADDITIONAL_LIBRARIES      "")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_ADDITIONAL_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}                         ${GCC_ADDITIONAL_LINK_FLAGS}" )

add_executable(puresu2gauge ${CODE_FILES_SU2})
add_executable(puresu3gauge ${CODE_FILES_SU3})

target_link_libraries(puresu2gauge ${GCC_ADDITIONAL_LIBRARIES})
target_link_libraries(puresu3gauge ${GCC_ADDITIONAL_LIBRARIES})

