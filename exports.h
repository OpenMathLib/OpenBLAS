#ifndef EXPORTS_H
#define EXPORTS_H

#pragma once

/***
 *  For GCC/Clang, always use -fvisibility=hidden. Then mark exported function
 *   implementations with OPENBLAS_EXPORT. For MSVC, using `__declspec` onces makes
 *    the default atrribute for any function in the entire shared object hidden
 *     (observed behaviour, documentation source needed).
 *     **/
#if defined (_WIN32) || defined (__CYGWIN__)
#  if defined (__GNUC__)
     /* GCC */
#    define OPENBLAS_EXPORT __attribute__ ((dllexport))
#    define OPENBLAS_IMPORT __attribute__ ((dllimport))
#  else
     /* MSVC */
#    define OPENBLAS_EXPORT __declspec(dllexport)
#    define OPENBLAS_IMPORT __declspec(dllimport)
#  endif
#else
   /* All other platforms. */
#  define OPENBLAS_EXPORT __attribute__ ((visibility ("default")))
#  define OPENBLAS_IMPORT
#endif


#endif // EXPORTS_H
