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

#include "ntt-twisted-radix-2-basic.hpp"
// #include "ntt-barrett.hpp"
#include "ntt-classical-radix-2-basic.hpp"

#include "ntt-twisted-radix-2-avx.hpp"
#include "ntt-classical-radix-2-avx.hpp"

template <static_modint_concept ModT>
void ntt_twisted(std::span<ModT> f) {
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::ntt_twisted_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::ntt_twisted_avx<ModT, false>(f);
    else
      detail::ntt_twisted_avx<ModT, true>(f);
  } else if constexpr (raw32_modint_concept<ModT>) {
    detail::ntt_twisted_basic(f);
  } else {
    detail::ntt_twisted_basic(f);
  }
}

template <static_modint_concept ModT>
void intt_twisted(std::span<ModT> f) {
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::intt_twisted_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::intt_twisted_avx<ModT, false>(f);
    else
      detail::intt_twisted_avx<ModT, true>(f);
  } else if constexpr (raw32_modint_concept<ModT>) {
    detail::intt_twisted_basic(f);
  } else {
    detail::intt_twisted_basic(f);
  }
}

template <static_modint_concept ModT>
void ntt_classical(std::span<ModT> f) {
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::ntt_classical_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::ntt_classical_avx<ModT, false>(f);
    else
      detail::ntt_classical_avx<ModT, true>(f);
  } else if constexpr (raw32_modint_concept<ModT>) {
    detail::ntt_classical_basic(f);
  } else {
    detail::ntt_classical_basic(f);
  }
}

template <static_modint_concept ModT>
void intt_classical(std::span<ModT> f) {
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::intt_classical_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::intt_classical_avx<ModT, false>(f);
    else
      detail::intt_classical_avx<ModT, true>(f);
  } else if constexpr (raw32_modint_concept<ModT>) {
    detail::intt_classical_basic(f);
  } else {
    detail::intt_classical_basic(f);
  }
}

template <static_modint_concept ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  ntt_classical(f);
  // ntt_twisted(f);
}

template <static_modint_concept ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  intt_classical(f);
  // intt_twisted(f);
}

#endif // ALGO_MATH_POLY_NTT
