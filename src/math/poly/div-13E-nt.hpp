#ifndef ALGO_MATH_POLY_DIV13ENT
#define ALGO_MATH_POLY_DIV13ENT

#include "poly-def.hpp"
#include "inv-10E-nt.hpp"

template <class ModT>
AVec<ModT> poly_div_13E(std::span<const ModT> lhs, std::span<const ModT> rhs, u32 m) {
  if (lhs.empty() || rhs.empty())
    return {};
  if (m == 1)
    return {lhs[0] / rhs[0]};
  u32 t = std::bit_ceil(m);
  AVec<ModT> g(t), x(t), q(t);
  AVec<ModT> u = poly_inv_10E<ModT>(rhs, t / 2); // 5E
  std::copy(rhs.begin(), std::min(rhs.begin() + t, rhs.end()), g.begin());
  std::copy(lhs.begin(), std::min(lhs.begin() + t / 2, lhs.end()), x.begin());
  u.resize(t);

  ntt<ModT>(x); // 1E
  ntt<ModT>(u); // 1E
  dot<ModT>(x, u);
  intt<ModT>(x); // 1E
  std::copy_n(x.begin(), t / 2, q.begin());
  ntt<ModT>(q); // 1E
  ntt<ModT>(g); // 1E
  dot<ModT>(q, g);
  intt<ModT>(q); // 1E
  std::fill_n(q.begin(), t / 2, 0);
  for (u32 i = t / 2; i < std::min<u32>(t, lhs.size()); ++i)
    q[i] = lhs[i] - q[i];
  ntt<ModT>(q); // 1E
  dot<ModT>(q, u);
  intt<ModT>(q); // 1E
  std::copy_n(q.begin() + t / 2, t / 2, x.begin() + t / 2);
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_DIV13ENT
