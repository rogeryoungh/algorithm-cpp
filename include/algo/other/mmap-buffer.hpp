#pragma once

#include "../base.hpp"

#include <cstdio>
#include <sys/mman.h>
#include <sys/stat.h>

ALGO_BEGIN_NAMESPACE

struct MmapBuf {
  struct stat sb;
  std::FILE *const f;
  const u8 *p, *beg, *end;
  MmapBuf(std::FILE *const f, usize) : f(f) {
    i32 fd = fileno(f);
    fstat(fd, &sb);
    beg = (u8 *)mmap(nullptr, sb.st_size + 4, PROT_READ, MAP_PRIVATE, fd, 0);
    p = beg, end = p + sb.st_size;
    madvise(const_cast<u8 *>(beg), sb.st_size + 4, MADV_SEQUENTIAL);
  }
  ~MmapBuf() {
    munmap(const_cast<u8 *>(beg), sb.st_size + 4);
  }
  bool eof() const {
    return end <= p;
  }
  void reserve(usize) {}
  u8 top() const {
    return *p;
  }
  u8 pop() {
    return *p++;
  }
};

ALGO_END_NAMESPACE
