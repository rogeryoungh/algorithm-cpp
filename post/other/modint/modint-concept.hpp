#ifndef ALGO_MODINT_CONCEPT
#define ALGO_MODINT_CONCEPT

#include <type_traits>

template <class ModT>
concept static_modint_concept = ModT::isStatic::value;

template <class ModT>
concept raw32_modint_concept = ModT::rawU32::value;

template <class ModT>
concept static_raw32_modint_concept = static_modint_concept<ModT> && raw32_modint_concept<ModT>;

template <class ModT>
concept runtime_modint_concept = !
ModT::isStatic::value;

template <class ModT>
concept montgomery_modint_concept = ModT::isMontgomery::value;

template <class ModT>
concept static_basic_modint_concept = !
montgomery_modint_concept<ModT> &&static_modint_concept<ModT>;

#endif // ALGO_MODINT_CONCEPT
