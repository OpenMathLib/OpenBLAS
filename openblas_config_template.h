/*This is only for "make install" target.*/

#ifdef NEEDBUNDERSCORE
#define BLASFUNC(FUNC) FUNC##_
#else
#define BLASFUNC(FUNC) FUNC
#endif

#ifdef QUAD_PRECISION
typedef struct {
  unsigned long x[2];
}  xdouble;
#elif defined EXPRECISION
#define xdouble long double
#else
#define xdouble double
#endif

#if defined(OS_WINDOWS) && defined(__64BIT__)
typedef long long BLASLONG;
typedef unsigned long long BLASULONG;
#else
typedef long BLASLONG;
typedef unsigned long BLASULONG;
#endif

#ifdef USE64BITINT
typedef BLASLONG blasint;
#else
typedef int blasint;
#endif

#if defined(XDOUBLE) || defined(DOUBLE)
#define FLOATRET	FLOAT
#else
#ifdef NEED_F2CCONV
#define FLOATRET	double
#else
#define FLOATRET	float
#endif
#endif
