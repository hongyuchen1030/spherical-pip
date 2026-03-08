#pragma once

#include <array>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <vector>

// Shewchuk wrapper header used by orient3d.hpp for SPIP_PREDICATES_REAL + exactinit().
#include "predicates.h"

// Config: defines FPG_UNCERTAIN_VALUE and pulls in fabs()
#include "spip/predicates/pck_filter_config.hpp"

// Vendored PCK filter (you placed it in third_party/quadruple3d.h).
// Your build already exposes third_party/ as an include dir.
#include "quadruple3d.h"

// Exact fallback implementation uses GEO::expansion_nt in a .cpp.
namespace spip::predicates {

enum class Sign : int8_t { Negative = -1, Zero = 0, Positive = 1 };

namespace internal {

using Real = SPIP_PREDICATES_REAL;

inline void ensure_exactinit() {
  static std::once_flag flag;
  std::call_once(flag, []() { exactinit(); });
}

template <class T>
inline constexpr bool supported_scalar_v =
    std::is_same_v<T, float> || std::is_same_v<T, double>;

template <class T>
inline void require_supported_and_matches_real() {
  static_assert(supported_scalar_v<T>, "spip::predicates::quadruple3d: only float/double are supported");
  static_assert(std::is_same_v<T, Real>,
                "spip::predicates::quadruple3d: scalar type must match SPIP_PREDICATES_REAL "
                "(define/un-define SPIP_PREDICATES_USE_FLOAT consistently)");
}

inline void require_nonnull(const void* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(std::string("spip::predicates::quadruple3d: null pointer for ") + name);
  }
}

inline void require_size3(std::size_t n, const char* name) {
  if (n != 3) {
    throw std::invalid_argument(std::string("spip::predicates::quadruple3d: ") + name + " must have size 3");
  }
}

inline Sign sign_from_int(int v) {
  if (v > 0) return Sign::Positive;
  if (v < 0) return Sign::Negative;
  return Sign::Zero;
}

// Exact fallback implemented in src/predicates/quadruple3d_exact.cpp
Sign quadruple3d_exact_fallback(const double* a, const double* b,
                                const double* c, const double* d);

inline Sign quadruple3d_raw_ptr(const Real* a,
                                const Real* b,
                                const Real* c,
                                const Real* d) {
  require_nonnull(a, "a");
  require_nonnull(b, "b");
  require_nonnull(c, "c");
  require_nonnull(d, "d");

  // Keep same init pattern as orient3d.hpp (harmless even if not strictly needed here).
  ensure_exactinit();

#ifdef SPIP_PREDICATES_USE_FLOAT
  static_assert(!std::is_same_v<Real, float>,
                "spip::predicates::quadruple3d: float mode not supported yet (PCK filter is double).");
#endif

  // Filter is double-based; Real should be double in your build.
  const int s = quadruple_3d_filter(
      reinterpret_cast<const double*>(a),
      reinterpret_cast<const double*>(b),
      reinterpret_cast<const double*>(c),
      reinterpret_cast<const double*>(d));

  if (s != FPG_UNCERTAIN_VALUE) {
    return sign_from_int(s);
  }

  // Exact fallback using GEO::expansion_nt
  return quadruple3d_exact_fallback(
      reinterpret_cast<const double*>(a),
      reinterpret_cast<const double*>(b),
      reinterpret_cast<const double*>(c),
      reinterpret_cast<const double*>(d));
}

}  // namespace internal

// Pointer API
template <class T>
inline Sign quadruple3d(const T* a, const T* b, const T* c, const T* d) {
  internal::require_supported_and_matches_real<T>();
  return internal::quadruple3d_raw_ptr(a, b, c, d);
}

// std::array API
template <class T>
inline Sign quadruple3d(const std::array<T, 3>& a,
                        const std::array<T, 3>& b,
                        const std::array<T, 3>& c,
                        const std::array<T, 3>& d) {
  internal::require_supported_and_matches_real<T>();
  return internal::quadruple3d_raw_ptr(a.data(), b.data(), c.data(), d.data());
}

// std::vector API
template <class T>
inline Sign quadruple3d(const std::vector<T>& a,
                        const std::vector<T>& b,
                        const std::vector<T>& c,
                        const std::vector<T>& d) {
  internal::require_supported_and_matches_real<T>();
  internal::require_size3(a.size(), "a");
  internal::require_size3(b.size(), "b");
  internal::require_size3(c.size(), "c");
  internal::require_size3(d.size(), "d");
  return internal::quadruple3d_raw_ptr(a.data(), b.data(), c.data(), d.data());
}

}  // namespace spip::predicates