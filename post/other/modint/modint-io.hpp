#ifdef ALGO_FASTIO

#ifdef ALGO_MODINT_BASIC_MODINT
template <u32 MOD>
inline FastI &operator>>(FastI &is, BasicModint<MOD> &m) {
  i64 x;
  is >> x;
  m = BasicModint<MOD>(x);
  return is;
}

template <u32 MOD>
inline FastO &operator<<(FastO &os, const BasicModint<MOD> &m) {
  return os << m.val();
}
#endif // ALGO_MODINT_BASIC_MODINT

#endif // ALGO_FASTIO
