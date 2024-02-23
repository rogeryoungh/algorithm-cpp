#pragma once

#include "../base.hpp"

#include "../number/montgomery.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class U>
inline constexpr U legendre(const Mont<U> &M, TI<U> a) {
  return M.qpow(a, (M.MOD - 1) / 2);
}

template <class U>
constexpr std::make_signed_t<U> cipolla(const Mont<U> &_M, TI<U> x) {
  if (x == 0)
    return 0;
  const Mont<U> M = _M;
  U n = M.trans(x);
  if (M.ncmp(legendre(M, n), M.ONE))
    return -1;
  if (M.MOD == 2)
    return 1;
  U a = 0, I = 0;
  for (; a < M.MOD; ++a) {
    I = M.sub(M.mul(a, a), n);
    if (M.cmp(legendre(M, I), M.neg(M.ONE))) {
      break;
    }
  }
  using Fp2 = std::array<U, 2>;
  Fp2 z = {M.ONE, M.ONE}, w = {a, M.ONE};
  auto fp2mul = [&](const Fp2 &x, const Fp2 &y) {
    U a = M.add(M.mul(x[0], y[0]), M.mul(M.mul(x[1], y[1]), I));
    U b = M.add(M.mul(x[1], y[0]), M.mul(y[1], x[0]));
    return Fp2{a, b};
  };
  for (U b = (M.MOD + 1) / 2; b; b /= 2) {
    if (b % 2 == 1)
      z = fp2mul(z, w);
    w = fp2mul(w, w);
  }
  U v1 = M.get(z[0]), v2 = M.MOD - v1;
  return v1 <= v2 ? v1 : v2;
}

ALGO_END_NAMESPACE
