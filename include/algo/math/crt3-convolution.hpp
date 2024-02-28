#pragma once

#include "../base.hpp"
#include "../number/montgomery.hpp"

#include <algorithm>
#include <span>
#include <vector>

ALGO_BEGIN_NAMESPACE

template <class NTT>
struct CRT3Convolution {
  u32 P1 = 167772161, G1 = 3;
  u32 P2 = 469762049, G2 = 3;
  u32 P3 = 754974721, G3 = 11;
  Mont32 M1{P1}, M2{P2}, M3{P3};
  NTT ntt1{M1, G1}, ntt2{M2, G2}, ntt3{M3, G3};
  u64 iv2 = M2.get(M2.inv(M2.trans(P1)));
  u64 iv3 = M3.get(M3.inv(M3.trans(P2)));
  u64 iv23 = M3.get(M3.inv(M3.mul(M3.trans(P1), M3.trans(P2))));
  std::vector<u32> conv(std::span<const u32> f, std::span<const u32> g, u32 mod) {
    u32 L = std::bit_ceil(f.size() + g.size() - 1);

    auto conv3 = [L, &f, &g](auto &&ntt) {
      std::vector<u32> f1(L), g1(L);
      std::ranges::copy(f, f1.begin());
      std::ranges::copy(g, g1.begin());
      ntt.conv(f1.data(), g1.data(), L);
      return f1;
    };

    auto f1 = conv3(ntt1);
    auto f2 = conv3(ntt2);
    auto f3 = conv3(ntt3);

    u64 M12 = u64(P1) * P2 % mod;
    for (u32 i = 0; i != L; ++i) {
      u64 a1 = f1[i], a2 = f2[i], a3 = f3[i];
      u64 x1 = (a2 - a1 + P2) * iv2 % P2;
      u64 x2 = ((a3 - a1 + P3) * iv23 + (P3 - x1) * iv3) % P3;
      f1[i] = (x2 * M12 + x1 * P1 + a1) % mod;
    }
    f1.resize(f.size() + g.size() - 1);
    return f1;
  }
};

ALGO_END_NAMESPACE
