#pragma once

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

struct Barrett64 {
  enum { s = 96 };
  static constexpr u128 s2 = u128(1) << s;
  u32 m;
  u128 ivm;
  Barrett64(u32 m_) : m(m_), ivm((s2 - 1) / m + 1) {}
  u32 div(u64 a) const {
    return a * ivm >> s;
  }
  u32 calc(u64 a) const {
    return a - u64(div(a)) * m;
  }
};

ALGO_END_NAMESPACE
