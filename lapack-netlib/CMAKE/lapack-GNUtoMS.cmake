# Skip conversion for non-GNU tools.
if(MINGW OR MSYS OR CYGWIN)
  return()
endif()

# Replace each imported target's import library.
foreach(lib ${ALL_TARGETS})
  # Replace for all imported build configurations.
  get_property(configs TARGET ${lib} PROPERTY IMPORTED_CONFIGURATIONS)
  foreach(config ${configs})
    get_property(implib TARGET ${lib} PROPERTY IMPORTED_IMPLIB_${config})
    # Switch to the MS-compatible import library.
    string(REGEX REPLACE "\\.dll\\.a$" ".lib" implib "${implib}")
    set_property(TARGET ${lib} PROPERTY IMPORTED_IMPLIB_${config} ${implib})
  endforeach()
endforeach()
