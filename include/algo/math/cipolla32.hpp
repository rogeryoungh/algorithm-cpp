#pragma once

#include "../base.hpp"

#include "../number/mont32.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

inline constexpr m32 legendre32(const Mont32 &M, m32 a) {
  return M.qpow(a, (M.MOD - 1) / 2);
}

constexpr i32 cipolla32(const Mont32 &_M, u32 x) {
  if (x == 0)
    return 0;
  const Mont32 M = _M;
  m32 n = M.trans(x);
  if (M.ncmp(legendre32(M, n), M.ONE))
    return -1;
  if (M.MOD == 2)
    return 1;
  u32 a = 0, I = 0;
  for (; a < M.MOD; ++a) {
    I = M.sub(M.mul(a, a), n);
    if (M.cmp(legendre32(M, I), M.neg(M.ONE))) {
      break;
    }
  }
  using Fp2 = std::array<m32, 2>;
  Fp2 z = {M.ONE, M.ONE}, w = {a, M.ONE};
  auto fp2mul = [&](const Fp2 &x, const Fp2 &y) {
    m32 a = M.add(M.mul(x[0], y[0]), M.mul(M.mul(x[1], y[1]), I));
    m32 b = M.add(M.mul(x[1], y[0]), M.mul(y[1], x[0]));
    return Fp2{a, b};
  };
  for (u32 b = (M.MOD + 1) / 2; b; b /= 2) {
    if (b % 2 == 1)
      z = fp2mul(z, w);
    w = fp2mul(w, w);
  }
  u32 v1 = M.get(z[0]), v2 = M.MOD - v1;
  return v1 <= v2 ? v1 : v2;
}

ALGO_END_NAMESPACE
