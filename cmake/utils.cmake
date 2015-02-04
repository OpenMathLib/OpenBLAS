# Functions to help with the OpenBLAS build

# Reads string from getarch into CMake vars. Format of getarch vars is VARNAME=VALUE
function(ParseGetArchVars GETARCH_IN)
  string(REGEX MATCHALL "[0-9_a-zA-Z]+=[0-9_a-zA-Z]+" GETARCH_RESULT_LIST "${GETARCH_IN}")
  foreach (GETARCH_LINE ${GETARCH_RESULT_LIST})
    # split the line into var and value, then assign the value to a CMake var
    string(REGEX MATCHALL "[0-9_a-zA-Z]+" SPLIT_VAR "${GETARCH_LINE}")
    list(GET SPLIT_VAR 0 VAR_NAME)
    list(GET SPLIT_VAR 1 VAR_VALUE)
    set(${VAR_NAME} ${VAR_VALUE} PARENT_SCOPE)
  endforeach ()
endfunction ()

# Returns all combinations of the input list, as a list with colon-separated combinations
# E.g. input of A B C returns A B C A:B A:C B:C
# N.B. The input is meant to be a list, and to past a list to a function in CMake you must quote it (e.g. AllCombinations("${LIST_VAR}")).
function(AllCombinations list_in)
  list(LENGTH list_in list_count)
  set(num_combos 1)
  # subtract 1 since we will iterate from 0 to num_combos
  math(EXPR num_combos "(${num_combos} << ${list_count}) - 1")
  set(LIST_OUT "")
  foreach (c RANGE 0 ${num_combos})
    set(current_combo "")
    # this is a little ridiculous just to iterate through a list w/ indices
    math(EXPR last_list_index "${list_count} - 1")
    foreach (list_index RANGE 0 ${last_list_index})
      math(EXPR bit "1 << ${list_index}")
      math(EXPR combo_has_bit "${c} & ${bit}")
      list(GET list_in ${list_index} list_elem)
      if (combo_has_bit)
        if (current_combo)
          set(current_combo "${current_combo}:${list_elem}")
        else ()
          set(current_combo ${list_elem})
        endif ()
      endif ()
    endforeach ()
    list(APPEND LIST_OUT ${current_combo})
  endforeach ()
  list(APPEND LIST_OUT " ") # Empty set is a valic combination, but CMake isn't appending the empty string for some reason, use a space
  set(LIST_OUT ${LIST_OUT} PARENT_SCOPE)
endfunction ()

# generates object files for each of the sources for each of the combinations of the preprocessor definitions passed in
# @param sources_in the source files to build from
# @param defines_in the preprocessor definitions that will be combined to create the object files
# @param all_defines_in (optional) preprocessor definitions that will be applied to all objects
function(GenerateObjects sources_in defines_in all_defines_in)
  AllCombinations("${defines_in}")
  set(define_combos ${LIST_OUT})
  set(OBJ_LIST_OUT "")
  foreach (source_file ${sources_in})
    foreach (def_combo ${define_combos})

      # replace colon separated list with semicolons, this turns it into a CMake list that we can use foreach with
      string(REPLACE ":" ";" def_combo ${def_combo})

      # build a unique variable name for this obj file by picking two letters from the defines (can't use one in this case)
      set(obj_name "")
      foreach (combo_elem ${def_combo})
        string(REGEX MATCH "^[A-Z][A-Z]" letter ${combo_elem})
        set(obj_name "${obj_name}${letter}")
      endforeach ()

      # parse file name
      string(REGEX MATCH "^[a-zA-Z_0-9]+" source_name ${source_file})
      string(TOUPPER ${source_name} source_name)

      # prepend the uppercased file name to the obj name
      set(obj_name "${source_name}_${obj_name}_OBJS")

      # now add the object and set the defines
      add_library(${obj_name} OBJECT ${source_file})
      set(cur_defines ${def_combo})
      if ("${cur_defines}" STREQUAL " ")
        set(cur_defines ${all_defines_in})
      else ()
        list(APPEND cur_defines ${all_defines_in})
      endif ()
      if (cur_defines AND NOT "${cur_defines}" STREQUAL " ") # using space as the empty set
        set_target_properties(${obj_name} PROPERTIES COMPILE_DEFINITIONS "${cur_defines}")
      endif ()
      list(APPEND OBJ_LIST_OUT ${obj_name})
    endforeach ()
  endforeach ()
  set(OBJ_LIST_OUT ${OBJ_LIST_OUT} PARENT_SCOPE)
endfunction ()

# generates object files for each of the sources, using the BLAS naming scheme to pass the funciton name as a preprocessor definition
# @param sources_in the source files to build from
# @param float_type_in the float type to define for this build (e.g. SINGLE/DOUBLE/etc)
# @param defines_in (optional) preprocessor definitions that will be applied to all objects
function(GenerateNamedObjects sources_in float_type_in defines_in)
  set(OBJ_LIST_OUT "")
  foreach (source_file ${sources_in})

    get_filename_component(source_name ${source_file} NAME_WE)

    string(SUBSTRING ${float_type_in} 0 1 float_char)
    string(TOLOWER ${float_char} float_char)

    # build a unique variable name for this obj file by picking two letters from the defines (can't use one in this case)
    set(obj_name "${float_char}${source_name}")

    # parse file name
    string(REGEX MATCH "^[a-zA-Z_0-9]+" source_name ${source_file})
    string(TOUPPER ${source_name} source_name)

    # now add the object and set the defines
    add_library(${obj_name} OBJECT ${source_file})
    set(obj_defines "ASMNAME=${FU}${obj_name};ASMFNAME=${FU}${obj_name}${BU};NAME=${obj_name}${BU};CNAME=${obj_name};CAR_NAME=\"${obj_name}${BU}\";CHAR_CNAME=\"${obj_name}\"")
    list(APPEND obj_defines ${defines_in})
    list(APPEND obj_defines ${float_type_in})
    set_target_properties(${obj_name} PROPERTIES COMPILE_DEFINITIONS "${obj_defines}")

    list(APPEND OBJ_LIST_OUT ${obj_name})

  endforeach ()
  set(OBJ_LIST_OUT ${OBJ_LIST_OUT} PARENT_SCOPE)
endfunction ()
