#ifndef ALGO_MATH_POLY_BASE
#define ALGO_MATH_POLY_BASE

#include "poly-def.hpp"

#include "inv-10E-nt-block.hpp"
#include "div-10E-nt-block.hpp"
#include "ln.hpp"
#include "exp-14E-nt-block.hpp"
#include "sqrt-8E-nt-block.hpp"
#include "invsqrt-12E-nt.hpp"
#include "pow.hpp"
#include "safe-sqrt.hpp"
#include "safe-pow.hpp"
#include "deriv.hpp"
#include "integr.hpp"

template <class ModT>
class Poly : public AVec<ModT> {
  using Vec = AVec<ModT>;

public:
  using Vec::empty;
  using Vec::resize;
  using Vec::size;
  using Vec::operator[];
  using Vec::begin;
  using Vec::cbegin;
  using Vec::cend;
  using Vec::end;

  static constexpr auto m_inv = poly_inv_10E_block<ModT>;
  static constexpr auto m_invsqrt = poly_invsqrt_12E<ModT>;
  static constexpr auto m_deriv = poly_deriv<ModT>;
  static constexpr auto m_integr = poly_integr<ModT>;
  static constexpr auto m_div = poly_div_10E_block<ModT, m_inv>;
  static constexpr auto m_ln = poly_ln<ModT, m_div>;
  static constexpr auto m_exp = poly_exp_14E_block<ModT, m_inv>;
  static constexpr auto m_sqrt = poly_sqrt_8E_block<ModT, m_inv>;
  static constexpr auto m_safe_sqrt = poly_safe_sqrt<ModT, m_sqrt>;
  static constexpr auto m_pow = poly_pow<ModT, m_ln, m_exp>;
  static constexpr auto m_safe_pow = poly_safe_pow<ModT, m_pow>;

  Poly() = default;

  Poly(u32 len) : Vec(len) {}

  Poly(const std::vector<u32> &v) : Vec(v.begin(), v.end()) {}
  Poly(const std::vector<ModT> &v) : Vec(v.begin(), v.end()) {}

  Poly(AVec<ModT> v) : Vec(std::move(v)) {}

  Poly(const ModT &v) : Vec({v}) {}

  Poly &operator+=(const Poly &rhs) {
    if (size() < rhs.size())
      resize(rhs.size());
    for (u32 i = 0; i < rhs.size(); ++i)
      (*this)[i] += rhs[i];
    return *this;
  }

  Poly &operator-=(const Poly &rhs) {
    if (size() < rhs.size())
      resize(rhs.size());
    for (u32 i = 0; i < rhs.size(); ++i)
      (*this)[i] -= rhs[i];
    return *this;
  }

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

  friend Poly operator+(const Poly &lhs, const Poly &rhs) {
    return Poly(lhs) += rhs;
  }

  friend Poly operator-(const Poly &lhs, const Poly &rhs) {
    return Poly(lhs) -= rhs;
  }

  Poly inv(u32 m) const {
    return m_inv(*this, m);
  }

  Poly deriv(u32 m) const {
    return m_deriv(*this, m);
  }

  Poly integr(u32 m, u32 C = 0) const {
    return m_integr(*this, m, C);
  }

  Poly div(const Poly &rhs, u32 m) {
    return m_div(*this, rhs, m);
  }

  Poly ln(u32 m) const {
    return m_ln(*this, m);
  }

  Poly exp(u32 m) const {
    return m_exp(*this, m);
  }

  Poly invsqrt(u32 m) const {
    return m_invsqrt(*this, m, this->front().sqrt().value().inv());
  }

  Poly sqrt(u32 m) const {
    return m_sqrt(*this, m, this->front().sqrt().value());
  }

  std::optional<Poly> sqrt_safe(u32 m) const {
    return m_safe_sqrt(*this, m);
  }

  Poly pow(u64 k, u32 m) const {
    return m_pow(*this, k, m);
  }

  Poly safe_pow(u64 k, u64 k2, u32 m) const {
    return m_safe_pow(*this, k, k2, m);
  }
};

#endif // ALGO_MATH_POLY_BASE
