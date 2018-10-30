#ifdef CBLAS

#include <stdarg.h>
#include "common.h"

void CNAME(blasint p, char *rout, char *form, ...)
{
   va_list args;

   va_start(args, form);

   if (p)
      fprintf(stderr, "Parameter %d to routine %s was incorrect\n", p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

