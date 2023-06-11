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

/////////////////////

#ifndef ALGO_DISABLE_NTT_CLASSICAL

#ifndef ALGO_DISABLE_NTT_RADIX_4

#include "ntt-classical-radix-4-basic.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2

#include "ntt-classical-radix-4-avx.hpp"

// classical-radix-4-avx2

template <class ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::ntt_classical_basic4(f);
    else if (u64(f.data()) & 0x1f)
      detail::ntt_classical_avx4<ModT, false>(f);
    else
      detail::ntt_classical_avx4<ModT, true>(f);
  } else {
    detail::ntt_classical_basic4(f);
  }
}

template <class ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::intt_classical_basic4(f);
    else if (u64(f.data()) & 0x1f)
      detail::intt_classical_avx4<ModT, false>(f);
    else
      detail::intt_classical_avx4<ModT, true>(f);
  } else {
    detail::intt_classical_basic4(f);
  }
}

#else // ALGO_DISABLE_SIMD_AVX2

#include "ntt-classical-radix-4-basic.hpp"

// classical-radix-4-basic

template <class ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  detail::ntt_classical_basic4(f);
}

template <class ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  detail::intt_classical_basic4(f);
}

#endif // ALGO_DISABLE_SIMD_AVX2

#else // ALGO_DISABLE_NTT_RADIX_4

#include "ntt-classical-radix-2-basic.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2

#include "ntt-classical-radix-2-avx.hpp"

// classical-radix-2-avx2

template <class ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::ntt_classical_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::ntt_classical_avx<ModT, false>(f);
    else
      detail::ntt_classical_avx<ModT, true>(f);
  } else {
    detail::ntt_classical_basic(f);
  }
}

template <class ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::intt_classical_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::intt_classical_avx<ModT, false>(f);
    else
      detail::intt_classical_avx<ModT, true>(f);
  } else {
    detail::intt_classical_basic(f);
  }
}

#else // ALGO_DISABLE_SIMD_AVX2

#include "ntt-classical-radix-2-basic.hpp"

// classical-radix-2-basic

template <class ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  detail::ntt_classical_basic(f);
}

template <class ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  detail::intt_classical_basic(f);
}

#endif // ALGO_DISABLE_SIMD_AVX2

#endif // ALGO_DISABLE_NTT_RADIX_4

#else // ALGO_DISABLE_NTT_CLASSICAL

#include "ntt-twisted-radix-2-basic.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2

#include "ntt-twisted-radix-2-avx.hpp"

// twisted-radix-2-avx2

template <class ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::ntt_twisted_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::ntt_twisted_avx<ModT, false>(f);
    else
      detail::ntt_twisted_avx<ModT, true>(f);
  } else {
    detail::ntt_twisted_basic(f);
  }
}

template <class ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  if constexpr (montgomery_modint_concept<ModT>) {
    if (f.size() < 16)
      detail::intt_twisted_basic(f);
    else if (u64(f.data()) & 0x1f)
      detail::intt_twisted_avx<ModT, false>(f);
    else
      detail::intt_twisted_avx<ModT, true>(f);
  } else {
    detail::intt_twisted_basic(f);
  }
}

#else // ALGO_DISABLE_SIMD_AVX2

// twisted-radix-2-basic

template <class ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  detail::ntt_twisted_basic(f);
}

template <class ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit<u32>(f.size()));
  detail::ntt_size += f.size();
  detail::intt_twisted_basic(f);
}

#endif // ALGO_DISABLE_SIMD_AVX2

#endif // ALGO_DISABLE_NTT_CLASSICAL

#endif // ALGO_MATH_POLY_NTT
