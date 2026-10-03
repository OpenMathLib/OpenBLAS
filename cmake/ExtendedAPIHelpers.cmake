include_guard(GLOBAL)

set(EXTENDED_API_GENERATOR
  "${CMAKE_CURRENT_LIST_DIR}/GenerateSuffixedSource.cmake" CACHE INTERNAL
  "Script that generates suffixed sources for extended APIs")

if(BUILD_INDEX64_EXT_API)
  add_custom_target(64bit_codegen ALL
    COMMENT "Generating 64-bit suffixed sources for extended API")
endif()

# Generate suffixed sources for an extended API. The generation happens at
# build time in the GenerateSuffixedSource.cmake script. Arguments:
#   SUFFIX <suffix>            -- suffix appended to symbols and file names
#                                 (required)
#   SYMBOL_ALLOWLIST <sym>...  -- only rename the listed symbols
#                                 (case-insensitive); by default every symbol
#                                 defined in the sources is renamed
#   NO_STRING_REPLACEMENTS     -- do not rename symbols inside string literals
function(generate_suffixed_sources target source_list generated_sources)
  set(options NO_STRING_REPLACEMENTS)
  set(oneValueArgs SUFFIX)
  set(multiValueArgs SYMBOL_ALLOWLIST)
  cmake_parse_arguments(PARSE_ARGV 3 extended_api "${options}" "${oneValueArgs}" "${multiValueArgs}")

  if(NOT DEFINED extended_api_SUFFIX)
    message(FATAL_ERROR "generate_suffixed_sources: SUFFIX is required")
  endif()

  get_filename_component(destination "${target}${extended_api_SUFFIX}_sources" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  get_property(_generated_suffixed_source_files GLOBAL PROPERTY EXTENDED_API_GENERATED_SOURCE_FILES)
  set(new_generated_source_files)
  set(generated_source_files)

  foreach(_source IN LISTS ${source_list})
    get_filename_component(source_abs "${_source}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    get_filename_component(source_name "${_source}" NAME_WLE)
    get_filename_component(source_ext "${_source}" EXT)
    set(output_file "${destination}/${source_name}${extended_api_SUFFIX}${source_ext}")

    set(_fortran_extensions ".f" ".F" ".f90" ".F90")
    if(NOT source_ext IN_LIST _fortran_extensions)
      message(WARNING "Skipping non-Fortran source '${_source}' for target '${target}'")
      continue()
    endif()

    # Make sure we only have one custom command generating a given output file
    if(NOT "${output_file}" IN_LIST _generated_suffixed_source_files)
      set(generator_args
        "-DINPUT_FILE=${source_abs}"
        "-DOUTPUT_FILE=${output_file}"
        "-DSUFFIX=${extended_api_SUFFIX}")
      if(extended_api_NO_STRING_REPLACEMENTS)
        list(APPEND generator_args "-DREPLACE_IN_STRINGS=OFF")
      endif()
      if(DEFINED extended_api_SYMBOL_ALLOWLIST)
        # Join with $<SEMICOLON> so that the allowlist stays one -D argument
        # on the command line but still reaches the script as a CMake list.
        # An allowlist change re-runs the generation (the command changes).
        list(JOIN extended_api_SYMBOL_ALLOWLIST "$<SEMICOLON>" _symbol_allowlist)
        list(APPEND generator_args "-DSYMBOL_ALLOWLIST=${_symbol_allowlist}")
      endif()

      add_custom_command(
        OUTPUT "${output_file}"
        COMMAND
          "${CMAKE_COMMAND}"
          ${generator_args}
          -P "${EXTENDED_API_GENERATOR}"
        DEPENDS
          "${source_abs}"
          "${EXTENDED_API_GENERATOR}"
        COMMENT "Generating ${extended_api_SUFFIX} extended API source for ${_source}"
        VERBATIM)

      list(APPEND new_generated_source_files "${output_file}")
    endif()

    list(APPEND generated_source_files "${output_file}")
  endforeach()

  # Make sure each generated source file is only part of one target to
  # avoid multiple targets trying to generate the same file
  if(new_generated_source_files)
    add_custom_target("${target}${extended_api_SUFFIX}_codegen" ALL
      DEPENDS ${new_generated_source_files}
      COMMENT "Generating ${extended_api_SUFFIX} suffixed sources for target ${target}")
    if(extended_api_SUFFIX STREQUAL "_64" AND TARGET 64bit_codegen)
      add_dependencies(64bit_codegen "${target}${extended_api_SUFFIX}_codegen")
    endif()

    set_property(GLOBAL APPEND PROPERTY EXTENDED_API_GENERATED_SOURCE_FILES ${new_generated_source_files})
  endif()

  set(${generated_sources} ${generated_source_files} PARENT_SCOPE)
endfunction()

# Generate 64-bit suffixed sources for the extended API. Kept as a thin
# wrapper around generate_suffixed_sources for the Index-64 extended API.
function(generate_64bit_suffixed_sources target source_list generated_sources)
  set(options NO_STRING_REPLACEMENTS)
  cmake_parse_arguments(PARSE_ARGV 3 extended_api "${options}" "" "")

  set(_forwarded_options)
  if(extended_api_NO_STRING_REPLACEMENTS)
    list(APPEND _forwarded_options NO_STRING_REPLACEMENTS)
  endif()

  generate_suffixed_sources("${target}" "${source_list}" _generated_source_files
    SUFFIX "_64" ${_forwarded_options})

  set(${generated_sources} ${_generated_source_files} PARENT_SCOPE)
endfunction()
