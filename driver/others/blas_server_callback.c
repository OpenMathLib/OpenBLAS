#include "common.h"

/* global variable to change threading backend from openblas-managed to caller-managed */
openblas_threads_callback openblas_threads_callback_ = 0;
void *openblas_threads_callback_data_ = 0;

/* non-threadsafe function should be called before any other
   openblas function to change how threads are managed */
void openblas_set_threads_callback(openblas_threads_callback callback, void *callback_data)
{
  openblas_threads_callback_ = callback;
  openblas_threads_callback_data_ = callback_data;
}