#ifdef ALGO_FASTIO

#ifdef ALGO_MODINT_STATIC_MODINT
template <class Space>
inline FastI &operator>>(FastI &is, StaticModint<Space> &m) {
  i64 x;
  is >> x;
  m = StaticModint<Space>(x);
  return is;
}

template <class Space>
inline FastO &operator<<(FastO &os, const StaticModint<Space> &m) {
  return os << m.val();
}
#endif // ALGO_MODINT_STATIC_MODINT

#endif // ALGO_FASTIO
