// Windows-specific code, include with:
// #ifdef _MSC_VER
// #include <msvc_compat.h>
// #endif

//#define _CRT_SECURE_NO_WARNINGS
#ifndef MSVC_COMPAT_H
#define MSVC_COMPAT_H

#define inline __inline
#define __func__ __FUNCTION__

#ifndef drand48
#include <drand48.h>
#endif

#include <float.h>
#define isnan _isnan
static int isinf(double x) {
    int y = _finite(x);
    if(y == 0) {
       return 0;
    } else {
       return 1;
    }
}
#define alloca _alloca
#define atoll _atoi64

#define ftello ftell

typedef int8_t bool;

static int R_OK = 4;
#endif
