#pragma once
// CUDA 13.1 declares rsqrt/rsqrtf without the noexcept added by glibc 2.43.
// Include the standard C math declarations before CUDA's implicit runtime
// preinclude, temporarily suppressing only the new optional C23/GNU math set.
// GNU/POSIX features stay enabled for the C++ threading and locale headers.
#if defined(__linux__)
#include <features.h>
#if defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 43)) &&          \
    defined(__CUDACC_VER_MAJOR__) && __CUDACC_VER_MAJOR__ == 13 && __CUDACC_VER_MINOR__ == 1
#pragma push_macro("__USE_GNU")
#pragma push_macro("__GLIBC_USE_ISOC23")
#undef __USE_GNU
#undef __GLIBC_USE_ISOC23
#define __GLIBC_USE_ISOC23 0
#include <math.h>
#pragma pop_macro("__GLIBC_USE_ISOC23")
#pragma pop_macro("__USE_GNU")
#endif
#endif
