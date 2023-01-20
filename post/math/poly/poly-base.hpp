#ifndef ALGO_MATH_POLY_BASE
#define ALGO_MATH_POLY_BASE

#include "../../other/modint/modint-concept.hpp"
#include "ntt.hpp"
#include "vec-dots.hpp"

#include <vector>

template <static_modint_concept ModT>
class Poly : public std::vector<ModT> {
  using Vec = typename std::vector<ModT>;

public:
  using Vec::resize;
  using Vec::size;
  using Vec::operator[];
  using Vec::begin;
  using Vec::cbegin;
  using Vec::cend;
  using Vec::end;

  Poly() = default;

  Poly(u32 len) : Vec(len) {
    self.resize(len);
  }

  Poly(const std::vector<u32> &v) : Vec(v.begin(), v.end()) {}

  Poly(const std::vector<ModT> &v) : Vec(v.begin(), v.end()) {}

  Poly &operator*=(const Poly &rhs) {
    if (self.empty() || rhs.empty()) {
      return self.resize(0), self;
    } else {
      u32 m = self.size() + rhs.size() - 1;
      u32 n = std::bit_ceil(m);
      self.resize(n);
      NttInfo<ModT>::ntt(self);
      if (self.data() == rhs.data()) {
        VecDots<ModT>::dot(self, self);
      } else {
        Vec vr(n);
        std::copy(rhs.cbegin(), rhs.cend(), vr.begin());
        NttInfo<ModT>::ntt(vr);
        VecDots<ModT>::dot(self, vr);
      }
      NttInfo<ModT>::intt(self);
      return self.resize(m), self;
    }
  }

  friend Poly operator*(const Poly &lhs, const Poly &rhs) {
    return Poly(lhs) *= rhs;
  }

  template <class U = u32>
  auto to_vec() const {
    return std::vector<U>{begin(), end()};
  }
};

#endif // ALGO_MATH_POLY_BASE
