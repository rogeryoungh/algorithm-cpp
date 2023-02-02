#ifndef ALGO_MATH_POLY_NTT
#define ALGO_MATH_POLY_NTT

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

u32 ntt_size = 0;

} // namespace detail

// #include "ntt-basic.hpp"

// #include "ntt-barrett.hpp"

#include "ntt-avx.hpp"

#endif // ALGO_MATH_POLY_NTT
