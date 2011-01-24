#if defined(__PGI) || defined(__PGIC__)
COMPILER_PGI
#endif

#if defined(__PATHSCALE__) || defined(__PATHCC__)
COMPILER_PATHSCALE
#endif

#if defined(__INTEL_COMPILER) || defined(__ICC) || defined(__ECC)
COMPILER_INTEL
#endif

#if defined(__OPENCC__)
COMPILER_OPEN64
#endif

#if defined(__SUNPRO_C)
COMPILER_SUN
#endif

#if defined(__IBMC__) || defined(__xlc__)
COMPILER_IBM
#endif

#if defined(__DECCC__)
COMPILER_DEC
#endif

#if defined(__GNUC__)
COMPILER_GNU
#endif

#if defined(__linux__)
OS_LINUX
#endif

#if defined(__FreeBSD__)
OS_FreeBSD
#endif

#if defined(__NetBSD__)
OS_NetBSD
#endif

#if defined(__sun)
OS_SunOS
#endif

#if defined(__APPLE__)
OS_Darwin
#endif

#if defined(_AIX)
OS_AIX
#endif

#if defined(__OSF)
OS_OSF
#endif

#if defined(__WIN32) || defined(__WIN64) || defined(__WINNT)
OS_WINNT
#endif

#if defined(__CYGWIN__)
OS_CYGWIN
#endif

#if defined(__INTERIX)
OS_INTERIX
#endif

#if defined(__i386) || defined(_X86)
ARCH_X86
#endif

#if defined(__x86_64__) || defined(__amd64__)
ARCH_X86_64
#endif

#if defined(__powerpc___) || defined(__PPC__) || defined(_POWER)
ARCH_POWER
#endif

#ifdef __mips64
ARCH_MIPS64
#endif

#if defined(__mips32) || defined(__mips)
ARCH_MIPS32
#endif

#ifdef __alpha
ARCH_ALPHA
#endif

#if defined(__sparc) || defined(__sparc__)
ARCH_SPARC
#endif

#if defined(__ia64__) || defined(__ia64)
ARCH_IA64
#endif

#if defined(__LP64) || defined(__LP64__) || defined(__ptr64) || defined(__x86_64__) || defined(__amd64__) || defined(__64BIT__)
BINARY_64
#endif
