#pragma once
// Internal support header for the hostlineage transcription (not shipped:
// src-internal, like mlj_engine.hpp). Types and table declarations for
// src/hostlineage.cpp (kernels) and src/hostlineage_tables.cpp (vendored
// glibc-2.39 table data). Provenance and licensing: see
// include/kalos/detail/hostlineage.hpp and the CI-FDLIBM section of
// docs/contracts/x13-automatic-pipeline-contract.md.

#include <cstdint>

namespace kalos::detail::hostlineage_src {

// The transcription is validated bit-exactly against the little-endian host
// libm (glibc-2.39, x86-64); the byte-swapped IBM-APML table branches are
// deliberately not vendored, so a big-endian build must fail closed.
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
              "hostlineage kernels are vendored little-endian only");
#else
#error "hostlineage requires a little-endian target (see hostlineage.hpp provenance)"
#endif

using int4 = int;  // glibc sysdeps/ieee754/dbl-64/mydefs.h

// mydefs.h word order for little-endian targets (the only ones this
// transcription serves; see the endianness guard above).
#define HIGH_HALF 1
#define LOW_HALF 0

// C++ list-initialization forbids silent narrowing of unsigned hex literals
// into int (glibc's C sources have no such restriction; bit patterns above
// INT_MAX are common in the tables). IW() performs the same two's-complement
// wraparound GCC performs for C initializers.
#define IW(hex) static_cast<int>(hex##U)

union mynumber {
  int4 i[2];
  double x;
  double d;
};

union sincostab_t {
  int4 i[880];
  double x[440];
};

union asncs_t {
  int4 i[5136];
  double x[2568];
};

// glibc-2.39 math_config.h struct layouts (x86-64 build: __FP_FAST_FMA not
// defined, so log_data keeps tab2; EXP/LOG/POW table bits = 7/7/7).
struct exp_data {
  double invln2N;
  double shift;
  double negln2hiN;
  double negln2loN;
  double poly[4];  // Last four coefficients.
  double exp2_shift;
  double exp2_poly[5];
  double invlog10_2N;
  double neglog10_2hiN;
  double neglog10_2loN;
  double exp10_poly[5];
  uint64_t tab[2 * (1 << 7)];
};

struct log_data {
  double ln2hi;
  double ln2lo;
  double poly[5];   // First coefficient is 1.
  double poly1[11];
  struct { double invc, logc; } tab[1 << 7];
  struct { double chi, clo; } tab2[1 << 7];
};

struct pow_log_data {
  double ln2hi;
  double ln2lo;
  double poly[7];  // First coefficient is 1.
  struct { double invc, pad, logc, logctail; } tab[1 << 7];
};

extern const sincostab_t sincostab;  // sincostab.c (renamed __sincostab)
extern const asncs_t asncs;          // asincos.tbl
extern const double inroot[128];     // root.tbl
extern const double powtwo[28];      // powtwo.tbl
extern const mynumber cij[241][7];   // uatan.tbl
extern const struct exp_data exp_data;       // e_exp_data.c (renamed)
extern const struct log_data log_data;       // e_log_data.c (renamed)
extern const struct pow_log_data pow_log_data;  // e_pow_log_data.c (renamed)

}  // namespace kalos::detail::hostlineage_src
