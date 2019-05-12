cmake_minimum_required( VERSION 3.6 )

project( hdk_adaptive_viscosity )

# CMAKE_PREFIX_PATH must contain the path to the toolkit/cmake subdirectory of
# the Houdini installation. See the "Compiling with CMake" section of the HDK
# documentation for more details, which describes several options for
# specifying this path.
list( APPEND CMAKE_PREFIX_PATH "$ENV{HFS}/toolkit/cmake" )

# Locate Houdini's libraries and header files.
# Registers an imported library target named 'Houdini'.
find_package( Houdini REQUIRED )


# Add the cmake/Modules directory for external dependency location scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

if (USEEIGEN)
    add_definitions(-DUSEEIGEN)
    
    find_package(EIGEN3 REQUIRED)
    if (EIGEN3_FOUND)
        include_directories(${EIGEN3_INCLUDE_DIR})
    endif()

    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

if (USESINGLEPRECISION)
    add_definitions(-DUSESINGLEPRECISION)
endif()

add_subdirectory(Source)
