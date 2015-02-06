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
# #param absent_codes codes to use when an element is absent from a combination. For example, if you have TRANS;UNIT;UPPER you may want the code to be NNL when nothing is present.
# @returns LIST_OUT a list of combinations
#          CODES_OUT a list of codes corresponding to each combination, with N meaning the item is not present, and the first letter of the list item meaning it is presen
function(AllCombinations list_in absent_codes_in)
  list(LENGTH list_in list_count)
  set(num_combos 1)
  # subtract 1 since we will iterate from 0 to num_combos
  math(EXPR num_combos "(${num_combos} << ${list_count}) - 1")
  set(LIST_OUT "")
  set(CODES_OUT "")
  foreach (c RANGE 0 ${num_combos})

    set(current_combo "")
    set(current_code "")

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
        string(SUBSTRING ${list_elem} 0 1 code_char)
      else ()
        list(GET absent_codes_in ${list_index} code_char)
      endif ()
      set(current_code "${current_code}${code_char}")
    endforeach ()

    if (current_combo STREQUAL "")
      list(APPEND LIST_OUT " ") # Empty set is a valid combination, but CMake isn't appending the empty string for some reason, use a space
    else ()
      list(APPEND LIST_OUT ${current_combo})
    endif ()
    list(APPEND CODES_OUT ${current_code})

  endforeach ()

  set(LIST_OUT ${LIST_OUT} PARENT_SCOPE)
  set(CODES_OUT ${CODES_OUT} PARENT_SCOPE)
endfunction ()

# generates object files for each of the sources, using the BLAS naming scheme to pass the funciton name as a preprocessor definition
# @param sources_in the source files to build from
# @param float_type_in the float type to define for this build (e.g. SINGLE/DOUBLE/etc)
# @param defines_in (optional) preprocessor definitions that will be applied to all objects
# @param name_in (optional) if this is set this name will be used instead of the filename. Use a * to indicate where the float character should go, if no star the character will be prepended.
#                           e.g. with DOUBLE set, "i*max" will generate the name "idmax", and "max" will be "dmax"
# @param replace_last_with replaces the last character in the filename with this string (e.g. symm_k should be symm_TU)
# @param append_with appends the filename with this string (e.g. trmm_R should be trmm_RTUU or some other combination of characters)
function(GenerateNamedObjects sources_in float_type_in)

  if (DEFINED ARGV2)
    set(defines_in ${ARGV2})
  endif ()

  if (DEFINED ARGV3)
    set(name_in ${ARGV3})
  endif ()

  if (DEFINED ARGV4)
    set(use_cblas ${ARGV4})
  else ()
    set(use_cblas 0)
  endif ()

  if (DEFINED ARGV5)
    set(replace_last_with ${ARGV5})
  endif ()

  if (DEFINED ARGV6)
    set(append_with ${ARGV6})
  endif ()

  set(OBJ_LIST_OUT "")
  foreach (source_file ${sources_in})

    if (NOT float_type_in STREQUAL "")
      string(SUBSTRING ${float_type_in} 0 1 float_char)
      string(TOLOWER ${float_char} float_char)
    endif ()

    if (NOT name_in)
      get_filename_component(source_name ${source_file} NAME_WE)
      set(obj_name "${float_char}${source_name}")
    else ()
      # replace * with float_char
      if (${name_in} MATCHES "\\*")
        string(REPLACE "*" ${float_char} obj_name ${name_in})
      else ()
        set(obj_name "${float_char}${name_in}")
      endif ()
    endif ()

    if (replace_last_with)
      string(REGEX REPLACE ".$" ${replace_last_with} obj_name ${obj_name})
    else ()
      set(obj_name "${obj_name}${append_with}")
    endif ()

    # now add the object and set the defines
    set(obj_defines ${defines_in})

    if (use_cblas)
      set(obj_name "cblas_${obj_name}")
      list(APPEND obj_defines "CBLAS")
    endif ()

    list(APPEND obj_defines "ASMNAME=${FU}${obj_name};ASMFNAME=${FU}${obj_name}${BU};NAME=${obj_name}${BU};CNAME=${obj_name};CAR_NAME=\"${obj_name}${BU}\";CHAR_CNAME=\"${obj_name}\"")
    list(APPEND obj_defines ${defines_in})
    list(APPEND obj_defines ${float_type_in})

    add_library(${obj_name} OBJECT ${source_file})
    set_target_properties(${obj_name} PROPERTIES COMPILE_DEFINITIONS "${obj_defines}")

    list(APPEND OBJ_LIST_OUT ${obj_name})

  endforeach ()
  set(OBJ_LIST_OUT ${OBJ_LIST_OUT} PARENT_SCOPE)
endfunction ()

# generates object files for each of the sources for each of the combinations of the preprocessor definitions passed in
# @param sources_in the source files to build from
# @param defines_in the preprocessor definitions that will be combined to create the object files
# @param float_type_in the float type to define for this build (e.g. SINGLE/DOUBLE/etc)
# @param all_defines_in (optional) preprocessor definitions that will be applied to all objects
# @param replace_scheme If 1, replace the "k" in the filename with the define combo letters. E.g. symm_k.c with TRANS and UNIT defined will be symm_TU.
#                  If 0, it will simply append the code, e.g. symm_L.c with TRANS and UNIT will be symm_LTU.
#                  If 2, it will append the code with an underscore, e.g. symm.c with TRANS and UNIT will be symm_TU.
#                  If 3, it will insert the code *around* the last character with an underscore, e.g. symm_L.c with TRANS and UNIT will be symm_TLU (required by BLAS level2 objects).
#                  If 4, it will insert the code before the last underscore. E.g. trtri_U_parallel with TRANS will be trtri_UT_parallel
# @param alternate_name replaces the source name as the object name (define codes are still appended)
function(GenerateCombinationObjects sources_in defines_in absent_codes_in float_type_in all_defines_in replace_scheme)

  if (DEFINED ARGV6)
    set(alternate_name ${ARGV6})
  endif ()

  AllCombinations("${defines_in}" "${absent_codes_in}")
  set(define_combos ${LIST_OUT})
  set(define_codes ${CODES_OUT})

  set(COMBO_OBJ_LIST_OUT "")
  list(LENGTH define_combos num_combos)
  math(EXPR num_combos "${num_combos} - 1")

  foreach (c RANGE 0 ${num_combos})

    list(GET define_combos ${c} define_combo)
    list(GET define_codes ${c} define_code)

    foreach (source_file ${sources_in})

      # replace colon separated list with semicolons, this turns it into a CMake list that we can use foreach with
      string(REPLACE ":" ";" define_combo ${define_combo})

      # now add the object and set the defines
      set(cur_defines ${define_combo})
      if ("${cur_defines}" STREQUAL " ")
        set(cur_defines ${all_defines_in})
      else ()
        list(APPEND cur_defines ${all_defines_in})
      endif ()

      set(replace_code "")
      set(append_code "")
      if (replace_scheme EQUAL 1)
        set(replace_code ${define_code})
      else ()
        if (replace_scheme EQUAL 2)
          set(append_code "_${define_code}")
        elseif (replace_scheme EQUAL 3)
          # first extract the last letter
          string(REGEX MATCH "[a-zA-Z]\\." last_letter ${source_file})
          string(SUBSTRING ${last_letter} 0 1 last_letter) # remove period from match
          # break the code up into the first letter and the remaining (should only be 2 anyway)
          string(SUBSTRING ${define_code} 0 1 define_code_first)
          string(SUBSTRING ${define_code} 1 -1 define_code_second)
          set(replace_code "${define_code_first}${last_letter}${define_code_second}")
        elseif (replace_scheme EQUAL 4)
          # insert code before the last underscore and pass that in as the alternate_name
          get_filename_component(alternate_name ${source_file} NAME_WE)
          set(extra_underscore "")
          # check if filename has two underscores, insert another if not (e.g. getrs_parallel needs to become getrs_U_parallel not getrsU_parallel)
          string(REGEX MATCH "_[a-zA-Z]+_" underscores ${alternate_name})
          string(LENGTH "${underscores}" underscores)
          if (underscores EQUAL 0)
            set(extra_underscore "_")
          endif ()
          string(REGEX REPLACE "(.+)(_[^_]+)$" "\\1${extra_underscore}${define_code}\\2" alternate_name ${alternate_name})
        else()
          set(append_code ${define_code}) # replace_scheme should be 0
        endif ()
      endif ()

      GenerateNamedObjects("${source_file}" "${float_type_in}" "${cur_defines}" "${alternate_name}" 0 "${replace_code}" "${append_code}")
      list(APPEND COMBO_OBJ_LIST_OUT "${OBJ_LIST_OUT}")
    endforeach ()
  endforeach ()

  set(COMBO_OBJ_LIST_OUT ${COMBO_OBJ_LIST_OUT} PARENT_SCOPE)
endfunction ()

