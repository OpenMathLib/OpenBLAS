#ifndef CBLAS_TEST_FMT_H
#define CBLAS_TEST_FMT_H

/* BIFS = blasint format specifier */
#ifdef USE64BITINT
#if defined(OS_WINDOWS) && defined(__64BIT__)
#define BIFS "%lld"
#else
#define BIFS "%ld"
#endif
#else
#define BIFS "%d"
#endif

#endif	// CBLAS_TEST_FMT_H
