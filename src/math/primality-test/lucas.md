---
author: "rogeryoungh"
title: "Lucas 素性测试"
date: "2023-01-21"
---

## 理论

设 $D$ 为无平方整数，记 $O$ 为二次域 $\mathbb{Q}(\sqrt{D})$ 上的 [代数整数环](https://www.bananaspace.org/wiki/代数整数)，有类比的定理：

【二次域中的 Fermat 小定理】设 $a \in O$，$p$ 为奇素数，且 $p \nmid d$，则

$$
a^p = \begin{cases} a, &\text{若 } \left(\frac{D}{p}\right) = 1 \\ \overline{a}, &\text{若 } \left(\frac{D}{p}\right) = -1\end{cases} \pmod p
$$

由于乘幂并不容易计算，引入一类特殊的二阶递推序列：设 $P, Q \in \mathbb{Z}$，特征方程 $x^2 - Px+ Q = 0$ 的两根分别为 $a, b$，定义 $P, Q$ 对应的 Lucas 序列为

$$
U_n = \frac{a^n - b^n}{a - b},\quad V_n = a^n + b^n
$$

> 我没理解乘幂计算的在哪里有困难，不能快速幂吗？猜测是那个除 $2$ 比较头疼。

其与二次域的关系是

$$
\left(\frac{P + \sqrt{D}}{2}\right)^n = \frac{V_n + U_n \sqrt{D}}{2}, \quad D = P^2 - 4Q
$$

首先我们限定 $\gcd(QD, p) = 1$，保证 $a,b,a-b$ 在模 $p$ 下可逆。再记 $\varepsilon = \left(\frac{D}{p}\right)$，不难根据 Fermat 小定理证明 $U_{p - \varepsilon} \equiv 0 \pmod p$。

## 算法描述

输入奇数 $N$，选取 $P, Q$ 作为 Lucas 序列参数，使得 $\varepsilon = \left(\frac{D}{N}\right) \neq 0$。

1. 若 $U_{N - \varepsilon} \equiv 0 \pmod N$，则输出 $N$ 可能为奇素数。
2. 若 $U_{N - \varepsilon} \not\equiv 0 \pmod N$，则输出 $N$ 为合数。

由于 $2$ 在奇数下的逆元总是存在的，我们可以使用矩阵快速幂计算上式。

对于 $P,Q$ 的选取，有两种建议：

1. Selfridge 建议：取 $D$ 为 $5, -7, 9, -11,13,\cdots$ 中使得 $\left(\frac{D}{N}\right) = -1$ 的第一个数，选择 $P = 1, Q = \frac{1-D}{2}$。
2. Baillie 建议：取 $D$ 为 $5, 9, 13, 17, 21, \cdots$ 中使得使得 $\left(\frac{D}{N}\right) = -1$ 的第一个数，选择 $P$ 为大于 $\sqrt{D}$ 的最小奇数，$Q = \frac{P^2 - D}{4}$。

论文中提到，取到合适的 $D$ 的期望次数为 $2$。若尝试次数过多，则说明 $N$ 很可能是平方数，需要检测下。

## 备注

网上找到的很多 Lucas 检测代码都是关于 $n+1$ 的因式分解，我要是能分解还判素数干嘛。。

似乎大家都建议固定 $D$，然后 $(p, q) \to (p + 2, p + q + 1)$ 这样，但是通不过 [LOJ143](https://loj.ac/p/143)，而且反例很多的样子。。。反正计算 Jacobi 的复杂度是 $O(\log^2 n)$ 的，不如换 $D$。

## 参考

- [素数判定 - mathmu](http://mathmu.github.io/MTCAS/doc/PrimeTest.html)
- [如何编程判断一个数是否是质数？ - Wiener ES的回答](https://www.zhihu.com/question/308322307/answer/574767625)
- [Primality Proving 3.2 n+1 tests and the Lucas-Lehmer test](https://primes.utm.edu/prove/prove3_2.html)

## 代码

{{<code file="./lucas.hpp" lang="cpp">}}
