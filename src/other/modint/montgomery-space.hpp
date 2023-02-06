#ifndef ALGO_MODINT_MONTGOMERY_SPACE
#define ALGO_MODINT_MONTGOMERY_SPACE

#include "../../base.hpp"
#include <type_traits>

template <class T, T MOD>
struct MontgomerySpace;

template <u32 MOD>
struct MontgomerySpace<u32, MOD> {
  static_assert(2 < MOD && MOD < u32(1) << 30, "mod must in [3, 2^30)");
  static_assert(MOD % 2 == 1, "mod must be odd");

  using ValueT = u32;
  using TransT = u32;
  using rawU32 = std::false_type;
  using isMontgomery = std::true_type;

  constexpr static u32 get_nr() {
    u32 x = 1;
    for (i32 i = 0; i < 5; ++i)
      x *= 2 - x * MOD;
    return x;
  }

  consteval static u32 mod() {
    return MOD;
  }

  enum : u32 {
    R = u32(u64(1) << 32 % MOD),
    IR = u32(-get_nr()),
    MOD2 = MOD * 2,
  };

  constexpr static TransT trans(ValueT x) {
    return (u64(x) << 32) % MOD;
  }

  constexpr static u32 reduce(u64 x) {
    return (x + u64(u32(x) * IR) * MOD) >> 32;
  }

  constexpr static u32 reduce_m(u32 n) {
    return n >> 31 ? n + MOD : n;
  }

  constexpr static u32 reduce_2m(u32 n) {
    return n >> 31 ? n + MOD2 : n;
  }

  constexpr static u32 add(u32 a, u32 b) {
    return reduce_2m(a + b - MOD2);
  }

  constexpr static u32 sub(u32 a, u32 b) {
    return reduce_2m(a - b);
  }

  constexpr static u32 mul(u32 a, u32 b) {
    return reduce(u64(a) * b);
  }

  constexpr static u32 safe(i64 x) {
    return reduce_m(x % MOD);
  }

  constexpr static ValueT val(TransT x) {
    return reduce_m(reduce(x) - MOD);
  }

  constexpr static u32 shift2(u32 x) {
    x = reduce(x);
    return (x & 1 ? x + MOD : x) >> 1;
  }
};

#endif // ALGO_MODINT_MONTGOMERY_SPACE
