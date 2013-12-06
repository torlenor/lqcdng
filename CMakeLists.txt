cmake_minimum_required (VERSION 2.8)
project (LQCDNG)

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.

# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
# add_library (neib neib.cpp)
# add_library (chargeconfig chargeconfig.cpp

include_directories (${LQCDNG_SOURCE_DIR})

# Make sure the linker can find the Hello library once it is built.
link_directories (${LQCDNG_BINARY_DIR})

set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)

# set(CMAKE_CXX_FLAGS "-O3 -Wall -pg")        ## Optimize
set(CMAKE_CXX_FLAGS "-O3 -Wall")        ## Optimize
# set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
# add_executable (calcalljack calcalljack.cpp calcint.cpp calcjack.cpp)
add_executable (lqcdng lqcdng.cc su3.cc init.cc helper.cc simulation.cc)

# Link the executable to the Hello library.
# target_link_libraries (calcalljack calcint calcjack) 
