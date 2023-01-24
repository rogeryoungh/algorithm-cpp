#ifndef ALGO_MATH_POLY_BASE
#define ALGO_MATH_POLY_BASE

#include "../../other/modint/modint-concept.hpp"
#include "ntt.hpp"
#include "vec-dots.hpp"
#include "inv-10E-nt.hpp"
#include "div-13E-nt.hpp"
#include "ln-13E-nt.hpp"
#include "deriv.hpp"
#include "integr.hpp"

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

  Poly(const std::vector<u32> &v) : Vec(v) {}

  Poly(const std::vector<ModT> &v) : Vec(v) {}

  Poly &operator*=(const Poly &rhs) {
    if (self.empty() || rhs.empty()) {
      return self.resize(0), self;
    } else {
      u32 m = self.size() + rhs.size() - 1;
      u32 n = std::bit_ceil(m);
      self.resize(n);
      ntt<ModT>(self);
      if (self.data() == rhs.data()) {
        dot<ModT>(self, self);
      } else {
        Vec vr(n);
        std::copy(rhs.cbegin(), rhs.cend(), vr.begin());
        ntt<ModT>(vr);
        dot<ModT>(self, vr);
      }
      intt<ModT>(self);
      return self.resize(m), self;
    }
  }

  friend Poly operator*(const Poly &lhs, const Poly &rhs) {
    return Poly(lhs) *= rhs;
  }

  Poly inv(u32 m) const {
    return poly_inv_10E<ModT>(self, m);
  }

  Poly deriv(u32 m) const {
    return poly_deriv<ModT>(self, m);
  }

  Poly integr(u32 m, u32 C = 0) const {
    return poly_integr<ModT>(self, m, C);
  }

  Poly div(const Poly &rhs, u32 m) {
    return poly_div_13E<ModT>(self, rhs, m);
  }

  Poly ln(u32 m) {
    return poly_ln_13E<ModT>(self, m);
  }

  template <class U = u32>
  auto to_vec() const {
    return std::vector<U>{begin(), end()};
  }
};

#endif // ALGO_MATH_POLY_BASE
