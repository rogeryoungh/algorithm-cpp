#ifndef ALGO_MODINT_BASIC_MOD_SPACE
#define ALGO_MODINT_BASIC_MOD_SPACE

#include "../../base.hpp"
#include <type_traits>

template <class T, T MOD>
struct BasicModSpace;

template <u32 MOD>
struct BasicModSpace<u32, MOD> {
  static_assert(2 < MOD && MOD < u32(1) << 31, "mod must in [3, 2^31)");

  using ValueT = u32;
  using TransT = u32;
  using rawU32 = std::true_type;
  using isMontgomery = std::false_type;

  enum : u32 {
    MOD2 = MOD * 2,
  };

  constexpr static u32 mod() {
    return MOD;
  }

  constexpr static TransT trans(ValueT x) {
    return x;
  }

  static ValueT val(TransT x) {
    return reduce_m(x);
  }

  static u32 reduce_m(ValueT n) {
    return n >> 31 ? n + MOD : n;
  }

  static u32 reduce_2m(u32 n) {
    return n >> 31 ? n + MOD2 : n;
  }

  static u32 add(u32 a, u32 b) {
    return reduce_m(a + b - MOD);
  }

  static u32 sub(u32 a, u32 b) {
    return reduce_m(a - b);
  }

  static u32 mul(u32 a, u32 b) {
    return u64(a) * b % MOD;
  }

  static u32 safe(i64 x) {
    return reduce_m(x % MOD);
  }

  static u32 shift2(u32 x) {
    return (x & 1 ? x + MOD : x) >> 1;
  }
};

#endif // ALGO_MODINT_BASIC_MOD_SPACE
