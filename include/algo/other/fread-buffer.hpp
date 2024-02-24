#pragma once

#include "../base.hpp"

#include <cassert>
#include <cstdio>
#include <vector>

ALGO_BEGIN_NAMESPACE

struct FReadBuf {
  std::FILE *const f;
  std::vector<u8> buf;
  const u8 *p, *end;
  FReadBuf(std::FILE *const f, usize size) : f(f), buf(size) {
    assert(size >= 0x100);
    p = buf.data(), end = p + size - 8;
    usize cnt = std::fread(buf.data(), 1, end - p, f);
    buf[cnt] = 0;
  }
  void load() {
    u8 *p2 = std::move(p, end, buf.data());
    usize cnt = std::fread(p2, 1, end - p2, f);
    p2[cnt] = 0, p = buf.data();
  }
  void reserve(usize n) {
    if (end - p < i64(n))
      load();
  }
  u8 top() const {
    return *p;
  }
  u8 pop() {
    return *p++;
  }
};

ALGO_END_NAMESPACE
