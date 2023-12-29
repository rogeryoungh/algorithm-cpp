#ifndef ALGO_H_ARCH_USEAVX2
#define ALGO_H_ARCH_USEAVX2

#include "../../base.hpp"

#pragma GCC target("avx2")
#include <immintrin.h>

ALGO_BEGIN_NAMESPACE

using i256 = __m256i;
using i256u = __m256i_u;
using i32x8 = __m256i;
using i64x4 = __m256i;
using u32x8 = __m256i;
using u64x4 = __m256i;

ALGO_END_NAMESPACE

#endif