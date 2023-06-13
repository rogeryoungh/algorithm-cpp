#if defined(ALGO_FASTIO) && defined(ALGO_MODINT_STATIC_MODINT)
#ifndef ALGO_FASTIO_WITH_STATIC_MODINT
#define ALGO_FASTIO_WITH_STATIC_MODINT

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

#endif // ALGO_FASTIO_WITH_STATIC_MODINT
#endif

#if defined(ALGO_FASTIO) && defined(ALGO_MODINT_DYNAMIC_MODINT)

#ifndef ALGO_FASTIO_WITH_DYNAMIC_MODINT
template <class Space>
inline FastI &operator>>(FastI &is, DynamicModint<Space> &m) {
  i64 x;
  is >> x;
  m = DynamicModint<Space>(x);
  return is;
}

template <class Space>
inline FastO &operator<<(FastO &os, const DynamicModint<Space> &m) {
  return os << m.val();
}
#endif // ALGO_FASTIO_WITH_DYNAMIC_MODINT
#endif
