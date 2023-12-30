#ifndef ALGO_H_ARCH_USEAVX512F
#define ALGO_H_ARCH_USEAVX512F

#include "../../base.hpp"

#pragma GCC target("avx512f")
#include <immintrin.h>

ALGO_BEGIN_NAMESPACE

using i512 = __m512i;
using i512u = __m512i_u;
using i32x16 = __m512i;
using i64x8 = __m512i;
using u32x16 = __m512i;
using u64x8 = __m512i;

ALGO_END_NAMESPACE

#endif
