#pragma once

#include "../base.hpp"

#include <cassert>
#include <cstdio>
#include <vector>

ALGO_BEGIN_NAMESPACE

struct FreadBuf {
  std::FILE *const f;
  std::vector<u8> buf;
  const u8 *p, *end;
  usize read_cnt = 0;
  FreadBuf(std::FILE *const f, usize size) : f(f), buf(size) {
    assert(size >= 0x100);
    p = buf.data(), end = p + size - 8;
    read_cnt = std::fread(buf.data(), 1, end - p, f);
    buf[read_cnt] = 0;
  }
  void load() {
    u8 *p2 = std::move(p, end, buf.data());
    read_cnt = std::fread(p2, 1, end - p2, f);
    p2[read_cnt] = 0, p = buf.data();
  }
  bool eof() const {
    return read_cnt == 0;
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
