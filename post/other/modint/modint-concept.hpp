#ifndef ALGO_MODINT_CONCEPT
#define ALGO_MODINT_CONCEPT

#include <type_traits>

template <class ModT>
concept static_modint_concept = ModT::is_static::value;

template <class ModT>
concept raw32_modint_concept = ModT::raw_u32::value;

template <class ModT>
concept static_raw32_modint_concept = static_modint_concept<ModT> && raw32_modint_concept<ModT>;

template <class ModT>
concept runtime_modint_concept = !ModT::is_static::value;

#endif // ALGO_MODINT_CONCEPT
