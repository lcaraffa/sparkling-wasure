project ( LASexamples )

include_directories (BEFORE ../inc ../../LASzip/src )

# definitions for VC++
add_definitions ( -DUSE_AS_DLL -D_CRT_SECURE_NO_WARNINGS )

find_library ( LAS_LIBRARY las HINTS ../src/VC12 )



add_executable ( lasexample_write_only lasexample_write_only.cpp  )
target_link_libraries( lasexample_write_only ${LAS_LIBRARY} )


add_executable ( lasexample_add_rgb lasexample_add_rgb.cpp  )
target_link_libraries( lasexample_add_rgb ${LAS_LIBRARY} )
