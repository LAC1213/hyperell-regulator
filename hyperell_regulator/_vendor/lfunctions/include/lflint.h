// Copyright Edgar Costa 2026
// See LICENSE file for license details.
//
// lflint.h — minimal RAII C++ wrappers over the FLINT C API.
//
// Two value types in namespace lfun:
//   ZZ   — wraps fmpz_t       (arbitrary-precision integer)
//   ZZX  — wraps fmpz_poly_t  (integer polynomial in T)
//
// Covers exactly the surface used by this repo's examples; not a
// general-purpose flintxx replacement.

#pragma once

#include <ostream>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

namespace lfun {

class ZZ {
  fmpz_t z;
public:
  ZZ()                       { fmpz_init(z); }
  ZZ(slong n)                { fmpz_init(z); fmpz_set_si(z, n); }
  ZZ(ulong n)                { fmpz_init(z); fmpz_set_ui(z, n); }
  ZZ(const ZZ& o)            { fmpz_init(z); fmpz_set(z, o.z); }
  // Move: init the target to 0 first, then swap. This leaves the moved-from
  // object in a valid (zero) state that fmpz_clear can safely handle.
  ZZ(ZZ&& o) noexcept        { fmpz_init(z); fmpz_swap(z, o.z); }
  ~ZZ()                      { fmpz_clear(z); }

  ZZ& operator=(const ZZ& o)     { fmpz_set(z, o.z); return *this; }
  ZZ& operator=(ZZ&& o) noexcept { fmpz_swap(z, o.z); return *this; }
  ZZ& operator=(slong n)         { fmpz_set_si(z, n); return *this; }
  ZZ& operator=(ulong n)         { fmpz_set_ui(z, n); return *this; }

  fmpz*       _fmpz()       { return z; }
  const fmpz* _fmpz() const { return z; }

  bool is_zero() const { return fmpz_is_zero(z); }
  bool is_one()  const { return fmpz_is_one(z); }
  void set_zero()      { fmpz_zero(z); }
  void set_one()       { fmpz_one(z); }

  // unary
  ZZ operator-() const { ZZ r; fmpz_neg(r.z, z); return r; }

  // compound assignment
  ZZ& operator+=(const ZZ& o) { fmpz_add  (z, z, o.z); return *this; }
  ZZ& operator*=(const ZZ& o) { fmpz_mul  (z, z, o.z); return *this; }
  ZZ& operator*=(slong n)     { fmpz_mul_si(z, z, n); return *this; }
  ZZ& operator*=(ulong n)     { fmpz_mul_ui(z, z, n); return *this; }

  // Non-negative remainder (0 <= r < |modulus|), so e.g. disc % ZZ(4) stays
  // in 0..3 even for a negative discriminant ((-5) % 3 == 1, not -2).
  ZZ operator%(const ZZ& m) const { ZZ r; fmpz_mod   (r.z, z, m.z); return r; }
  ZZ operator%(ulong n)     const { ZZ r; fmpz_mod_ui(r.z, z, n);   return r; }

  bool operator==(const ZZ& o) const { return fmpz_equal   (z, o.z) != 0; }
  bool operator==(slong n)     const { return fmpz_equal_si(z, n)   != 0; }
  bool operator!=(const ZZ& o) const { return !(*this == o); }
  bool operator!=(slong n)     const { return !(*this == n); }
  bool operator< (const ZZ& o) const { return fmpz_cmp(z, o.z) < 0; }

  // fmpz_divisible_si is FLINT's signed-divisor variant; sign of n is ignored
  // for the divisibility test. divexact is exact division (UB if not divisible,
  // same contract as fmpz_divexact_si).
  bool divisible(slong n) const { return fmpz_divisible_si(z, n) != 0; }
  ZZ   divexact (slong n) const { ZZ r; fmpz_divexact_si(r.z, z, n); return r; }

  // Conversion. Supported T: slong, ulong, double. Any other T fires a
  // dependent static_assert when instantiated (loud compile error, not
  // a silent wrong value).
  template<class T> T to() const {
    static_assert(!sizeof(T*), "ZZ::to<T>: only slong, ulong, double are supported");
    return T{};  // unreachable; satisfies return-type checker
  }

  friend ZZ operator+(const ZZ&, const ZZ&);
  friend ZZ operator-(const ZZ&, const ZZ&);
  friend ZZ operator*(const ZZ&, const ZZ&);
  friend ZZ operator*(slong,     const ZZ&);
  friend ZZ operator*(const ZZ&, slong);
  friend ZZ pow(const ZZ&, ulong);
  friend std::ostream& operator<<(std::ostream&, const ZZ&);
};

inline ZZ operator+(const ZZ& a, const ZZ& b) { ZZ r; fmpz_add(r.z, a.z, b.z); return r; }
inline ZZ operator-(const ZZ& a, const ZZ& b) { ZZ r; fmpz_sub(r.z, a.z, b.z); return r; }
inline ZZ operator*(const ZZ& a, const ZZ& b) { ZZ r; fmpz_mul(r.z, a.z, b.z); return r; }
inline ZZ operator*(slong a,     const ZZ& b) { ZZ r; fmpz_mul_si(r.z, b.z, a); return r; }
inline ZZ operator*(const ZZ& a, slong b)     { ZZ r; fmpz_mul_si(r.z, a.z, b); return r; }

inline ZZ pow(const ZZ& a, ulong e) { ZZ r; fmpz_pow_ui(r.z, a.z, e); return r; }

inline std::ostream& operator<<(std::ostream& s, const ZZ& x) {
  char* str = fmpz_get_str(nullptr, 10, x.z);
  s << str;
  flint_free(str);
  return s;
}

template<> inline slong  ZZ::to<slong>()  const { return fmpz_get_si(z); }
template<> inline ulong  ZZ::to<ulong>()  const { return fmpz_get_ui(z); }
template<> inline double ZZ::to<double>() const { return fmpz_get_d (z); }

class ZZX {
  fmpz_poly_t p;
public:
  ZZX()                          { fmpz_poly_init(p); }
  explicit ZZX(slong alloc)      { fmpz_poly_init2(p, alloc); }
  ZZX(const ZZX& o)              { fmpz_poly_init(p); fmpz_poly_set(p, o.p); }
  ZZX(ZZX&& o) noexcept          { fmpz_poly_init(p); fmpz_poly_swap(p, o.p); }
  ~ZZX()                         { fmpz_poly_clear(p); }

  ZZX& operator=(const ZZX& o)     { fmpz_poly_set(p, o.p); return *this; }
  ZZX& operator=(ZZX&& o) noexcept { fmpz_poly_swap(p, o.p); return *this; }
  ZZX& operator=(slong n)          { fmpz_poly_set_si(p, n); return *this; }

  fmpz_poly_struct*       _poly()       { return p; }
  const fmpz_poly_struct* _poly() const { return p; }

  slong degree() const           { return fmpz_poly_degree(p); }
  void  fit_length(slong n)      { fmpz_poly_fit_length(p, n); }

  ZZ get_coeff(slong i) const {
    ZZ r;
    fmpz_poly_get_coeff_fmpz(r._fmpz(), p, i);
    return r;
  }
  void set_coeff(slong i, slong c)        { fmpz_poly_set_coeff_si(p, i, c); }
  void set_coeff(slong i, const ZZ& c)    { fmpz_poly_set_coeff_fmpz(p, i, c._fmpz()); }

  ZZX& operator+=(const ZZX& o) { fmpz_poly_add(p, p, o.p); return *this; }
  ZZX& operator*=(const ZZX& o) { fmpz_poly_mul(p, p, o.p); return *this; }
};

} // namespace lfun
