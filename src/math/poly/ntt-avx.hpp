#ifndef ALGO_MATH_POLY_NTT_AVX
#define ALGO_MATH_POLY_NTT_AVX

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

extern u32 ntt_size;

} // namespace detail

#include "../../other/avx.hpp"



template <static_modint_concept ModT>
auto &prepare_root_avx(u32 m) {
  using ValueT = typename ModT::ValueT;
  static constexpr ValueT P = ModT::mod();
  static constexpr ValueT g = 3;
  static constexpr ValueT max_bit = ValueT(1) << std::countr_zero(ModT::mod() - 1);

  static std::vector<ModT> rt{1, 1};
  assert(m <= max_bit);
  while (rt.size() < m) {
    u32 n = rt.size();
    rt.resize(n * 2);
    ModT p = ModT(g).pow((P - 1) / n / 2);
    for (u32 i = n; i < n * 2; i += 2) {
      rt[i] = rt[i / 2], rt[i + 1] = p * rt[i];
    }
  }
  return rt;
}

namespace detail {

template <montgomery_modint_concept ModT>
static void ntt_less_16(std::span<ModT> f) { // dif
  i32 n = f.size();
  auto &rt = prepare_root_avx<ModT>(n);

  for (i32 l = n / 2; l > 0; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = f[i + j + l];
        f[i + j + l] = rt[j + l] * (x - y);
        f[i + j] = x + y;
      }
    }
  }
}

template <montgomery_modint_concept ModT>
static void ntt_avx(std::span<ModT> f) { // dif
  i32 n = f.size();
  auto &rt = prepare_root_avx<ModT>(n);

  for (i32 l = n / 2; l >= 16; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; j += 8) {
        auto px = (u8x32 *)&f[i + j];
        auto py = (u8x32 *)&f[i + j + l];
        auto pr = (u8x32 *)&rt[j + l];
        m32x8 fx = i32x8(px), fy = i32x8(py);
        m32x8 r = i32x8(pr);

        m32x8 ry = r * (fx - fy);
        m32x8 rx = fx + fy;
        rx.v.store(px);
        ry.v.store(py);
      }
    }
  }
  auto gen_rt = [&](i32 L) {
    std::array<ModT, 8> ar;
    for (i32 i = 0; i < 8; i += L)
      for (i32 j = 0; j < L; ++j)
        ar[i + j] = rt[L + j];
    m32x8 r = i32x8(ar);
    return r;
  };
  auto r8 = gen_rt(8);
  auto r4 = gen_rt(4);
  auto r2 = gen_rt(2);
  // L = 8
  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];
    m32x8 fx = i32x8(px), fy = i32x8(py);

    m32x8 ry = r8 * (fx - fy);
    m32x8 rx = fx + fy;
    rx.v.store(px);
    ry.v.store(py);
  }
  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];

    // 01230123 45674567
    i32x8 fxy0 = px, fxy1 = py;
    // 01234567
    m32x8 fx = i32x8::permute_128<0x31>(fxy0, fxy1);
    m32x8 fy = i32x8::permute_128<0x20>(fxy0, fxy1);

    m32x8 ry = r4 * (fx - fy);
    m32x8 rx = fx + fy;
    i32x8::permute_128<0x20>(rx.v, ry.v).store(px);
    i32x8::permute_128<0x31>(rx.v, ry.v).store(py);
  }
  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];
    // 01012323 45456767
    i32x8 fxy0 = px, fxy1 = py;
    // 01452367
    m32x8 fx = i32x8::unpack_hi_64(fxy0, fxy1);
    m32x8 fy = i32x8::unpack_lo_64(fxy0, fxy1);

    m32x8 ry = r2 * (fx - fy);
    m32x8 rx = fx + fy;

    i32x8::unpack_hi_64(rx.v, ry.v).store(py);
    i32x8::unpack_lo_64(rx.v, ry.v).store(px);
  }
  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];
    // 00112233 44556677
    i32x8 fxy0 = px, fxy1 = py;
    fxy0 = fxy0.shuffle<0b11011000>();
    // 01012323 45456767
    fxy1 = fxy1.shuffle<0b11011000>();
    // 01452367
    m32x8 fx = i32x8::unpack_hi_64(fxy0, fxy1);
    m32x8 fy = i32x8::unpack_lo_64(fxy0, fxy1);

    m32x8 ry = fx - fy;
    m32x8 rx = fx + fy;

    i32x8::unpack_hi_64(rx.v, ry.v).shuffle<0b11011000>().store(py);
    i32x8::unpack_lo_64(rx.v, ry.v).shuffle<0b11011000>().store(px);
  }
}

template <montgomery_modint_concept ModT>
static void intt_less_16(std::span<ModT> f) { // dit
  i32 n = f.size();
  auto &rt = prepare_root_avx<ModT>(n);
  for (i32 l = 1; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = rt[j + l] * f[i + j + l];
        f[i + j] = x + y, f[i + j + l] = x - y;
      }
    }
  }
  const ModT ivn = ModT(n).inv();
  for (i32 i = 0; i < n; i++)
    f[i] *= ivn;
  std::reverse(f.begin() + 1, f.end());
}

template <montgomery_modint_concept ModT>
static void intt_avx(std::span<ModT> f) { // dit
  i32 n = f.size();
  auto &rt = prepare_root_avx<ModT>(n);
  auto gen_rt = [&](i32 L) {
    std::array<ModT, 8> ar;
    for (i32 i = 0; i < 8; i += L)
      for (i32 j = 0; j < L; ++j)
        ar[i + j] = rt[L + j];
    m32x8 r = i32x8(ar);
    return r;
  };
  auto r8 = gen_rt(8);
  auto r4 = gen_rt(4);
  auto r2 = gen_rt(2);
  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];
    // 00112233 44556677
    i32x8 fxy0 = px, fxy1 = py;
    fxy0 = fxy0.shuffle<0b11011000>();
    // 01012323 45456767
    fxy1 = fxy1.shuffle<0b11011000>();
    // 01452367
    m32x8 fy = i32x8::unpack_hi_64(fxy0, fxy1);
    m32x8 fx = i32x8::unpack_lo_64(fxy0, fxy1);

    m32x8 ry = fx - fy;
    m32x8 rx = fx + fy;

    i32x8::unpack_hi_64(rx.v, ry.v).shuffle<0b11011000>().store(py);
    i32x8::unpack_lo_64(rx.v, ry.v).shuffle<0b11011000>().store(px);
  }

  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];
    // 01012323 45456767
    i32x8 fxy0 = px, fxy1 = py;
    // 01452367
    m32x8 fy = i32x8::unpack_hi_64(fxy0, fxy1);
    m32x8 fx = i32x8::unpack_lo_64(fxy0, fxy1);

    fy *= r2;
    m32x8 ry = (fx - fy);
    m32x8 rx = fx + fy;

    i32x8::unpack_hi_64(rx.v, ry.v).store(py);
    i32x8::unpack_lo_64(rx.v, ry.v).store(px);
  }
  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];

    // 01230123 45674567
    i32x8 fxy0 = px, fxy1 = py;
    // 01234567
    m32x8 fx = i32x8::permute_128<0x20>(fxy0, fxy1);
    m32x8 fy = i32x8::permute_128<0x31>(fxy0, fxy1);

    fy *= r4;
    m32x8 ry = fx - fy;
    m32x8 rx = fx + fy;
    i32x8::permute_128<0x20>(rx.v, ry.v).store(px);
    i32x8::permute_128<0x31>(rx.v, ry.v).store(py);
  }
  for (i32 i = 0; i < n; i += 16) {
    auto px = (u8x32 *)&f[i];
    auto py = (u8x32 *)&f[i + 8];
    m32x8 fx = i32x8(px), fy = i32x8(py);

    fy *= r8;
    m32x8 ry = fx - fy;
    m32x8 rx = fx + fy;
    rx.v.store(px);
    ry.v.store(py);
  }
  for (i64 l = 16; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; j += 8) {
        auto px = (u8x32 *)&f[i + j];
        auto py = (u8x32 *)&f[i + j + l];
        auto pr = (u8x32 *)&rt[j + l];
        m32x8 r = i32x8(pr);
        m32x8 fx = i32x8(px), fy = m32x8(py) * r;
        m32x8 rx = fx + fy;
        m32x8 ry = fx - fy;
        rx.v.store(px);
        ry.v.store(py);
      }
    }
  }
  ModT ivn = ModT(n).inv();
  m32x8 ivn8 = i32x8(ivn.raw());
  for (i32 i = 0; i < n; i += 8) {
    auto pf = (u8x32 *)&f[i];
    m32x8 fi = i32x8((u8x32 *)&f[i]);
    fi *= ivn8;
    fi.v.store(pf);
  }
  std::reverse(f.begin() + 1, f.end());
}

} // namespace detail

template <montgomery_modint_concept ModT>
static void ntt(std::span<ModT> f) { // dif
  i32 n = f.size();
  assert(std::has_single_bit<u32>(n));
  detail::ntt_size += n;

  if (n < 16) {
    detail::ntt_less_16(f);
  } else {
    detail::ntt_avx(f);
  }
}

template <montgomery_modint_concept ModT>
static void intt(std::span<ModT> f) { // dit
  i32 n = f.size();
  assert(std::has_single_bit<u32>(n));
  detail::ntt_size += n;
  if (n < 16) {
    detail::intt_less_16(f);
  } else {
    detail::intt_avx(f);
  }
}

#endif // ALGO_MATH_POLY_NTT_AVX
