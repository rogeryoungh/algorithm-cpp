#ifndef ALGO_MODINT_DYNAMIC_MONTGOMERY_SPACE
#define ALGO_MODINT_DYNAMIC_MONTGOMERY_SPACE

#include "../../base.hpp"
#include <type_traits>

template <class T, i32 tag>
struct DynamicMontgomerySpace;

template <i32 tag>
struct DynamicMontgomerySpace<u32, tag> {

  using ValueT = u32;
  using TransT = u32;
  using rawU32 = std::false_type;
  using isMontgomery = std::true_type;

  inline static u32 R = -1;
  inline static u32 IR = -1;
  inline static u32 MOD = -1;
  inline static u32 MOD2 = -1;

  static bool set_mod(u32 mod) {
    if (mod <= 1 || mod >= (1u << 30))
      return false;
    if (mod % 2 == 0)
      return false;
    MOD = mod, MOD2 = MOD * 2;
    R = u64(1) << 32 % MOD;
    IR = -get_nr();
    return true;
  }

  static u32 get_nr() {
    u32 x = 1;
    for (i32 i = 0; i < 5; ++i)
      x *= 2 - x * MOD;
    return x;
  }

  static u32 mod() {
    return MOD;
  }

  static TransT trans(ValueT x) {
    return (u64(x) << 32) % MOD;
  }

  static u32 reduce(u64 x) {
    return (x + u64(u32(x) * IR) * MOD) >> 32;
  }

  static u32 reduce_m(u32 n) {
    return n >> 31 ? n + MOD : n;
  }

  static u32 reduce_2m(u32 n) {
    return n >> 31 ? n + MOD2 : n;
  }

  static u32 add(u32 a, u32 b) {
    return reduce_2m(a + b - MOD2);
  }

  static u32 sub(u32 a, u32 b) {
    return reduce_2m(a - b);
  }

  static u32 mul(u32 a, u32 b) {
    return reduce(u64(a) * b);
  }

  static u32 safe(i64 x) {
    return reduce_m(x % MOD);
  }

  static ValueT val(TransT x) {
    return reduce_m(reduce(x) - MOD);
  }

  static u32 shift2(u32 x) {
    x = reduce(x);
    return (x & 1 ? x + MOD : x) >> 1;
  }
};

#endif // ALGO_MODINT_MONTGOMERY_SPACE
