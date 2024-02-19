#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_mod"

#define ALGO_NO_NAMESPACE
#include "../../include/algo/math/cipolla32.hpp"

#include <iostream>

i32 main() {
  std::cin.tie(nullptr)->sync_with_stdio(false);

  u32 T;
  std::cin >> T;
  auto M9 = Mont32{998244353};
  while (T--) {
    u32 n, p;
    std::cin >> n >> p;
    if (p == 998244353) {
      std::cout << cipolla32(M9, n) << '\n';
    } else {
      std::cout << cipolla32(Mont32{p}, n) << '\n';
    }
  }
  return 0;
}
