#define PROBLEM "https://judge.yosupo.jp/problem/primality_test"

#define ALGO_NO_NAMESPACE
#include "../../include/algo/math/miller-rabin.hpp"

#include <iostream>

i32 main() {
  std::cin.tie(nullptr)->sync_with_stdio(false);

  u32 Q;
  std::cin >> Q;
  while (Q--) {
    u64 n;
    std::cin >> n;
    std::cout << (miller_rabin(n) ? "Yes\n" : "No\n");
  }
  return 0;
}
