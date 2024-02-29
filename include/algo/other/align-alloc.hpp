#pragma once

#include "../base.hpp"

#include <cstdlib>
#include <vector>

#ifndef ALGO_AVEC_ALIGN
#define ALGO_AVEC_ALIGN 32
#endif

ALGO_BEGIN_NAMESPACE

template <class T, u32 align>
struct AlignedAllocator {
  using value_type = T;
  T *allocate(std::size_t n) const {
    std::size_t size = sizeof(T) * n;
    void *ptr;
    if (size < align) {
      ptr = std::malloc(size);
    } else {
      size = (size + align - 1) & (~(align - 1));
      ptr = std::aligned_alloc(align, size);
    }
    return reinterpret_cast<T *>(ptr);
  }
  template <class U>
  struct rebind {
    using other = AlignedAllocator<U, align>;
  };
  void deallocate(T *p, std::size_t) const {
    std::free(p);
  }
};

template <class T, u32 align = ALGO_AVEC_ALIGN>
using AVec = std::vector<T, AlignedAllocator<T, align>>;

ALGO_END_NAMESPACE
