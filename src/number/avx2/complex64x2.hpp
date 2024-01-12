#ifndef ALGO_H_MATH_COMPLEX64X2
#define ALGO_H_MATH_COMPLEX64X2

#include "../../base.hpp"
#include "../../other/avx2/use-avx2.hpp"
#include "../complex64.hpp"

ALGO_BEGIN_NAMESPACE

struct CP64x2 {
  f64x4 v;
  CP64x2() : v() {}

  CP64x2(const f64x4 &o) : v(o) {}

  operator f64x4() const {
    return v;
  }

  static f64x4 from(CP64 v) {
    return _mm256_broadcast_pd(reinterpret_cast<__m128d *>(&v));
  }

  CP64x2 &operator+=(const CP64x2 &o) {
    v = _mm256_add_pd(v, o);
    return *this;
  }
  CP64x2 &operator-=(const CP64x2 &o) {
    v = _mm256_sub_pd(v, o);
    return *this;
  }
  // https://www.researchgate.net/figure/Vectorized-complex-multiplication-using-AVX-2_fig2_337532904
  CP64x2 &operator*=(const CP64x2 &o) {
    f64x4 cc = _mm256_shuffle_pd(o, o, 0x00);
    f64x4 dd = _mm256_shuffle_pd(o, o, 0x0f);
    f64x4 ba = _mm256_shuffle_pd(v, v, 0x05);
    v = _mm256_fmaddsub_pd(v, cc, _mm256_mul_pd(ba, dd));
    return *this;
  }
  static CP64x2 cmul(const CP64x2 &l, const CP64x2 &r) {
    f64x4 cc = _mm256_shuffle_pd(r, r, 0x00);
    f64x4 dd = _mm256_shuffle_pd(r, r, 0x0f);
    f64x4 ba = _mm256_shuffle_pd(l, l, 0x05);
    return _mm256_fmsubadd_pd(l, cc, _mm256_mul_pd(ba, dd));
  }
  CP64x2 &operator*=(f64 o) {
    v = _mm256_mul_pd(v, _mm256_set1_pd(o));
    return *this;
  }
  CP64x2 &operator/=(f64 o) {
    v = _mm256_div_pd(v, _mm256_set1_pd(o));
    return *this;
  }

  friend CP64x2 operator+(const CP64x2 &l, const CP64x2 &r) {
    return CP64x2(l) += r;
  }
  friend CP64x2 operator-(const CP64x2 &l, const CP64x2 &r) {
    return CP64x2(l) -= r;
  }
  friend CP64x2 operator*(const CP64x2 &l, const CP64x2 &r) {
    return CP64x2(l) *= r;
  }
  friend CP64x2 operator*(const CP64x2 &l, const f64 &r) {
    return CP64x2(l) *= r;
  }
  friend CP64x2 operator/(const CP64x2 &l, const f64 &r) {
    return CP64x2(l) /= r;
  }
};

ALGO_END_NAMESPACE

#endif
