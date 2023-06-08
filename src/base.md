---
author: "rogeryoungh"
title: "算法库"
date: "2023-01-20"
---

实现的算法列表。

## 数学相关

### 基础数论

- [Miller Rabin 素性测试](../math/primality-test/miller-rabin)，`math/primality-test/miller-rabin.hpp`。
- [Lucas 素性测试](../math/primality-test/lucas)，`math/primality-test/lucas.hpp`。
- [Baillie PSW 素性测试](../math/primality-test/baillie-psw)，`math/primality-test/baillie-psw.hpp`。

### Poly

- [Poly 骨架](../math/poly/poly)，`math/poly/poly-base.hpp`。
- [NTT](../math/poly/ntt)，`math/poly/ntt.hpp`。

#### 牛顿迭代

- [多项式乘法逆 10E](../math/poly/inv-10e-nt)，`math/poly/inv-10E-nt.hpp`。
- [多项式乘法逆 12E](../math/poly/inv-12e-nt)，`math/poly/inv-12E-nt.hpp`。
- [多项式除法 13E](../math/poly/div-13e-nt)，`math/poly/div-13E-nt.hpp`。
- [多项式 EXP 17E](../math/poly/exp-17e-nt)，`math/poly/exp-17E-nt.hpp`。
- [多项式开根 11E](../math/poly/sqrt-11e-nt)，`math/poly/sqrt-11E-nt.hpp`。

#### 分块牛顿迭代

- [多项式乘法逆 10E 分块](../math/poly/inv-10e-nt-block)，`math/poly/inv-10E-block.hpp`。
- [多项式除法 10E 分块](../math/poly/div-10e-nt-block)，`math/poly/div-10E-nt-block.hpp`。
- [多项式开根 8E 分块](../math/poly/sqrt-8e-nt-block)，`math/poly/sqrt-8E-nt-block.hpp`。
- [多项式 EXP 14E 分块](../math/poly/exp-14e-nt-block)，`math/poly/exp-17E-nt-block.hpp`。


## 数据结构相关

- [跳表](../datastruct/skiplist)`datastruct/skiplist.hpp`