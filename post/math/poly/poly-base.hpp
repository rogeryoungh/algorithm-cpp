#ifndef ALGO_MATH_POLY_BASE
#define ALGO_MATH_POLY_BASE

#include "../../other/modint/modint-concept.hpp"
#include "ntt.hpp"
#include "vec-dots.hpp"
#include "inv-10E-nt-block.hpp"
#include "div-10E-nt-block.hpp"
#include "ln.hpp"
#include "exp-14E-nt-block.hpp"
#include "sqrt-8E-nt-block.hpp"
#include "safe-sqrt.hpp"
#include "deriv.hpp"
#include "integr.hpp"

#include <optional>
#include <vector>

template <static_modint_concept ModT>
class Poly : public std::vector<ModT> {
  using Vec = typename std::vector<ModT>;

public:
  using Vec::empty;
  using Vec::resize;
  using Vec::size;
  using Vec::operator[];
  using Vec::begin;
  using Vec::cbegin;
  using Vec::cend;
  using Vec::end;

  Poly() = default;

  Poly(u32 len) : Vec(len) {}

  Poly(const std::vector<u32> &v) : Vec(v) {}

  Poly(const std::vector<ModT> &v) : Vec(v) {}

  Poly &operator*=(const Poly &rhs) {
    if (empty() || rhs.empty()) {
      return resize(0), *this;
    } else {
      u32 m = size() + rhs.size() - 1;
      u32 n = std::bit_ceil(m);
      resize(n);
      ntt<ModT>(*this);
      if (this->data() == rhs.data()) {
        dot<ModT>(*this, *this);
      } else {
        Vec vr(n);
        std::copy(rhs.cbegin(), rhs.cend(), vr.begin());
        ntt<ModT>(vr);
        dot<ModT>(*this, vr);
      }
      intt<ModT>(*this);
      return resize(m), *this;
    }
  }

  friend Poly operator*(const Poly &lhs, const Poly &rhs) {
    return Poly(lhs) *= rhs;
  }

  Poly inv(u32 m) const {
    return poly_inv_10E<ModT>(*this, m);
  }

  Poly deriv(u32 m) const {
    return poly_deriv<ModT>(*this, m);
  }

  Poly integr(u32 m, u32 C = 0) const {
    return poly_integr<ModT>(*this, m, C);
  }

  Poly div(const Poly &rhs, u32 m) {
    return poly_div_10E_block<ModT>(*this, rhs, m);
  }

  Poly ln(u32 m) const {
    return poly_ln<ModT>(*this, m, poly_div_10E_block<ModT>);
  }

  Poly exp(u32 m) const {
    return poly_exp_14E_block<ModT>(*this, m);
  }

  Poly sqrt(u32 m) const {
    return poly_sqrt_8E_block<ModT>(*this, m, this->front().sqrt());
  }

  std::optional<Poly> sqrt_safe(u32 m) const {
    return poly_safe_sqrt<ModT>(*this, m, poly_sqrt_8E_block<ModT>);
  }

  template <class U = u32>
  auto to_vec() const {
    return std::vector<U>{begin(), end()};
  }
};

#endif // ALGO_MATH_POLY_BASE
