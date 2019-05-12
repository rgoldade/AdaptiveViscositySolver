set( library_name HDK_AdaptiveViscosity )

# Add a library and its source files.
add_library( ${library_name} SHARED
    HDK_AdaptiveViscosity.cpp
    HDK_AdaptiveViscosity.h
    HDK_OctreeGrid.cpp
    HDK_OctreeGrid.h
    HDK_OctreeVectorFieldInterpolator.cpp
    HDK_OctreeVectorFieldInterpolator.h
)

# Link against the Houdini libraries, and add required include directories and
# compile definitions.
target_link_libraries( ${library_name} Houdini )

# Include ${CMAKE_CURRENT_BINARY_DIR} for the generated header.
target_include_directories( ${library_name} PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Sets several common target properties, such as the library's output directory.
houdini_configure_target( ${library_name} )
