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

#include "ntt-basic.hpp"

// #include "ntt-barrett.hpp"

#include "ntt-avx.hpp"

template <static_modint_concept ModT>
void ntt(std::span<ModT> f) {
  i32 n = f.size();
  assert(std::has_single_bit<u32>(n));
  detail::ntt_size += n;

  if constexpr (montgomery_modint_concept<ModT>) {
    if (n < 16)
      ntt_basic(f);
    else if (u64(f.data()) & 0x1f)
      ntt_avx<ModT, false>(f);
    else
      ntt_avx<ModT, true>(f);
  } else if constexpr (raw32_modint_concept<ModT>) {
    ntt_basic(f);
  } else {
    ntt_basic(f);
  }
}

template <static_modint_concept ModT>
void intt(std::span<ModT> f) {
  i32 n = f.size();
  assert(std::has_single_bit<u32>(n));
  detail::ntt_size += n;

  if constexpr (montgomery_modint_concept<ModT>) {
    if (n < 16)
      intt_basic(f);
    else if (u64(f.data()) & 0x1f)
      intt_avx<ModT, false>(f);
    else
      intt_avx<ModT, true>(f);
  } else if constexpr (raw32_modint_concept<ModT>) {
    intt_basic(f);
  } else {
    intt_basic(f);
  }
}

#endif // ALGO_MATH_POLY_NTT
