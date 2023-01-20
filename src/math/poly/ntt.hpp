#ifndef ALGO_MATH_POLY_NTT
#define ALGO_MATH_POLY_NTT

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

template <static_modint_concept ModT>
class NttInfo {
protected:
  enum : typename ModT::value_type {
    P = ModT::get_mod(),
    g = 3, // 原根
    max_bit = (P - 1) & -(P - 1),
    iv2 = (P + 1) / 2
  };
  inline static std::vector<ModT> rt{1, 1};
  inline static std::array<ModT, 32> ivn{1, iv2};

  static ModT get_ivn(u32 n) {
    return ivn[std::countr_zero(n)];
  }

  static void prepare_root(u32 m) {
    assert(m <= max_bit);
    while (rt.size() < m) {
      u32 n = rt.size();
      rt.resize(n * 2);
      ModT p = ModT(g).pow((P - 1) / n / 2);
      for (u32 i = n; i < n * 2; i += 2) {
        rt[i] = rt[i / 2], rt[i + 1] = p * rt[i];
      }
      ivn[std::countr_zero(n * 2)] = P - (P - 1) / n / 2;
    }
  }

public:
  static void ntt(std::span<ModT> f) { // dif
    u32 n = f.size();
    assert(std::has_single_bit(n));
    prepare_root(n), ntt_size += n;
    for (u32 l = n / 2; l > 0; l /= 2) {
      for (u32 i = 0; i < n; i += l * 2) {
        for (u32 j = 0; j < l; ++j) {
          ModT x = f[i + j], y = f[i + j + l];
          f[i + j + l] = (x - y) * rt[j + l];
          f[i + j] = x + y;
        }
      }
    }
  }

  static void intt(std::span<ModT> f) { // dit
    u32 n = f.size();
    assert(std::has_single_bit(n));
    prepare_root(n), ntt_size += n;
    f[n] = 1;
    for (u32 l = 1; l < n; l *= 2) {
      for (u32 i = 0; i < n; i += l * 2) {
        for (u32 j = 0; j < l; ++j) {
          ModT x = f[i + j], y = rt[j + l] * f[i + j + l];
          f[i + j] = x + y, f[i + j + l] = x - y;
        }
      }
    }
    const ModT ivn = get_ivn(n);
    for (u32 i = 0; i < n; i++)
      f[i] *= ivn;
    std::reverse(f.begin() + 1, f.end());
  }
};

#endif // ALGO_MATH_POLY_NTT
