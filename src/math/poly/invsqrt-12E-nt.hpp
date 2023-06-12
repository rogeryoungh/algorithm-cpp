#ifndef ALGO_MATH_POLY_INVSQRT_12ENT
#define ALGO_MATH_POLY_INVSQRT_12ENT

#include "poly-def.hpp"

template <class ModT>
auto poly_invsqrt_12E(std::span<const ModT> self, u32 m, const ModT &x0) {
  u32 n = std::bit_ceil(m);
  AVec<ModT> x(n * 2);
  x[0] = x0;
  ModT ivn2 = -ModT(2).inv();
  for (u32 t = 1; t < n; t *= 2) {
    AVec<ModT> u(t * 4), s(t * 4);
    std::copy(self.begin(), std::min(self.begin() + t * 2, self.end()), u.begin());
    std::copy(x.begin(), x.begin() + t, s.begin());
    ntt<ModT>(u); // 4E
    ntt<ModT>(s); // 4E
#ifndef ALGO_DISABLE_SIMD_AVX2
    auto ivn2x8 = simd::M32x8<ModT>::from(ivn2);
    vectorization_2<ModT, true>(t * 4, u.data(), s.data(), [ivn2, ivn2x8]<class T>(T &ui, T si) {
      if constexpr (std::is_same_v<T, ModT>)
        ui *= si * si * si * ivn2;
      else
        ui *= si * si * si * ivn2x8;
    });
#else
    vectorization_2<ModT, true>(t * 4, u.data(), s.data(), [ivn2]<class T>(T &ui, T si) {
      ui *= si * si * si * ivn2;
    });
#endif
    intt<ModT>(u); // 4E
    std::copy_n(u.begin() + t, t, x.begin() + t);
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INVSQRT_12ENT
