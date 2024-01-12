#ifndef ALGO_H_MATH_COMPLEX64
#define ALGO_H_MATH_COMPLEX64

#include "../base.hpp"
#include <cmath>

ALGO_BEGIN_NAMESPACE

struct CP64 {
  f64 x, y;
  explicit constexpr CP64(f64 x = 0, f64 y = 0) : x(x), y(y) {}
  constexpr CP64 &operator+=(const CP64 &o) {
    x += o.x, y += o.y;
    return *this;
  }
  constexpr CP64 &operator-=(const CP64 &o) {
    x -= o.x, y -= o.y;
    return *this;
  }
  constexpr CP64 &operator*=(const CP64 &o) {
    f64 r = x * o.x - y * o.y, i = x * o.y + y * o.x;
    x = r, y = i;
    return *this;
  }
  constexpr CP64 &operator*=(f64 o) {
    x *= o, y *= o;
    return *this;
  }
  constexpr CP64 &operator/=(f64 o) {
    x /= o, y /= o;
    return *this;
  }

  constexpr f64 norm() const {
    return x * x + y * y;
  }

  constexpr CP64 conj() const {
    return CP64{x, -y};
  }

  constexpr f64 arg() const {
    return std::atan2(y, x);
  }

  constexpr CP64 mulj() const {
    return CP64{-y, x};
  }

  constexpr CP64 &operator/=(const CP64 &o) {
    f64 r = x * o.x + y * o.y, i = -y * o.x + x * o.y;
    f64 n = o.norm();
    x = r / n, y = i / n;
    return *this;
  }

  constexpr friend CP64 operator+(const CP64 &l, const CP64 &r) {
    return CP64(l) += r;
  }
  constexpr friend CP64 operator-(const CP64 &l, const CP64 &r) {
    return CP64(l) -= r;
  }
  constexpr friend CP64 operator*(const CP64 &l, const CP64 &r) {
    return CP64(l) *= r;
  }
  constexpr friend CP64 operator/(const CP64 &l, const CP64 &r) {
    return CP64(l) /= r;
  }
  constexpr friend CP64 operator*(const CP64 &l, const f64 &r) {
    return CP64(l) *= r;
  }
  constexpr friend CP64 operator/(const CP64 &l, const f64 &r) {
    return CP64(l) /= r;
  }
  constexpr static CP64 cmul(const CP64 &l, const CP64 &r) {
    f64 x = l.x * r.x + l.y * r.y, y = l.y * r.x - l.x * r.y;
    return CP64{x, y};
  }

  constexpr static CP64 polar(f64 t) {
    return CP64{std::cos(t), std::sin(t)};
  }

  constexpr static CP64 polar(f64 t, f64 r) {
    return CP64{r * std::cos(t), r * std::sin(t)};
  }
};

ALGO_END_NAMESPACE

#endif
