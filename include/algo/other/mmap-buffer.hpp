#pragma once

#include "../base.hpp"

#include <cstdio>
#include <sys/mman.h>
#include <sys/stat.h>

ALGO_BEGIN_NAMESPACE

struct MmapBuf {
  struct stat sb;
  std::FILE *const f;
  u8 *p, *p0;
  MmapBuf(std::FILE *const f, usize) : f(f) {
    i32 fd = fileno(f);
    fstat(fd, &sb);
    p0 = p = (u8 *)mmap(nullptr, sb.st_size + 4, PROT_READ, MAP_PRIVATE, fd, 0);
    madvise(p, sb.st_size, MADV_SEQUENTIAL);
  }
  ~MmapBuf() {
    munmap(p0, sb.st_size + 4);
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
