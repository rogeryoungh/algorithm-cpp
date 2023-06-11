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

// classical-radix-4
#define ALGO_DETAIL_NTT detail::ntt_classical_basic4
#define ALGO_DETAIL_INTT detail::intt_classical_basic4
#include "ntt-classical-radix-4-basic.hpp"
#ifndef ALGO_DISABLE_SIMD_AVX2
#define ALGO_DETAIL_NTT_AVX detail::ntt_classical_avx4
#define ALGO_DETAIL_INTT_AVX detail::intt_classical_avx4
#include "ntt-classical-radix-4-avx.hpp"
#endif // ALGO_DISABLE_SIMD_AVX2

#else

// classical-radix-2
#define ALGO_DETAIL_NTT detail::ntt_classical_basic
#define ALGO_DETAIL_INTT detail::intt_classical_basic
#include "ntt-classical-radix-2-basic.hpp"
#ifndef ALGO_DISABLE_SIMD_AVX2
#define ALGO_DETAIL_NTT_AVX detail::ntt_classical_avx
#define ALGO_DETAIL_INTT_AVX detail::intt_classical_avx
#include "ntt-classical-radix-2-avx.hpp"
#endif // ALGO_DISABLE_SIMD_AVX2

#endif // ALGO_DISABLE_NTT_RADIX_4

#else

// twisted-radix-2
#define ALGO_DETAIL_NTT detail::ntt_twisted_basic
#define ALGO_DETAIL_INTT detail::intt_twisted_basic
#include "ntt-twisted-radix-2-basic.hpp"
#ifndef ALGO_DISABLE_SIMD_AVX2
#define ALGO_DETAIL_NTT_AVX detail::ntt_twisted_avx
#define ALGO_DETAIL_INTT_AVX detail::intt_twisted_avx
#include "ntt-twisted-radix-2-avx.hpp"
#endif // ALGO_DISABLE_SIMD_AVX2

#endif

template <class ModT>
void ntt(std::span<ModT> f) {
  assert(std::has_single_bit(f.size()));
  detail::ntt_size += f.size();
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && f.size() > 16) {
    if (u64(f.data()) & 0x1f)
      ALGO_DETAIL_NTT_AVX<ModT, false>(f);
    else
      ALGO_DETAIL_NTT_AVX<ModT, true>(f);
  } else {
#endif
    ALGO_DETAIL_NTT(f);
#ifndef ALGO_DISABLE_SIMD_AVX2
  }
#endif
}

template <class ModT>
void intt(std::span<ModT> f) {
  assert(std::has_single_bit(f.size()));
  detail::ntt_size += f.size();
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && f.size() > 16) {
    if (u64(f.data()) & 0x1f)
      ALGO_DETAIL_INTT_AVX<ModT, false>(f);
    else
      ALGO_DETAIL_INTT_AVX<ModT, true>(f);
  } else {
#endif
    ALGO_DETAIL_INTT(f);
#ifndef ALGO_DISABLE_SIMD_AVX2
  }
#endif
}

#undef ALGO_DETAIL_NTT
#undef ALGO_DETAIL_NTT_AVX
#undef ALGO_DETAIL_INTT
#undef ALGO_DETAIL_INTT_AVX

#endif // ALGO_MATH_POLY_NTT
