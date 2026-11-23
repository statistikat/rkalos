// Vendored glibc-2.39 dbl-64 transcendental kernels (host-lineage, CI-FDLIBM
// LEG A). Full provenance, pinned digests, lineage decisions, the exhaustive
// list of source transformations, and the licensing flag live in
// include/kalos/detail/hostlineage.hpp — read it before editing this file.
//
// The kernels below are transcribed from:
//   sysdeps/ieee754/dbl-64/s_sin.c + branred.c   (sin/cos, IBM APML)
//   sysdeps/ieee754/dbl-64/e_exp.c               (exp, ARM table-driven)
//   sysdeps/ieee754/dbl-64/e_log.c               (log, ARM table-driven)
//   sysdeps/ieee754/dbl-64/e_pow.c               (pow, ARM table-driven)
//   sysdeps/ieee754/dbl-64/e_hypot.c             (hypot, Borges MyHypot3)
//   sysdeps/ieee754/dbl-64/e_atan2.c             (atan2, IBM APML uatan2)
//   sysdeps/ieee754/dbl-64/e_asin.c              (acos/asin, IBM APML)
// Table data: src/hostlineage_tables.cpp (generated; see generator docstring).
//
// Fidelity rules observed here: every arithmetic expression keeps glibc's
// exact operation order and literal constants; the target compiles with
// -ffp-contract=off so no FMA contraction can occur; the only textual
// changes are (1) glibc plumbing macros replaced by the bit-equivalent
// helpers below (documented in the public header), (2) statement-expression
// macros converted to verbatim inline functions, (3) the namespace scoping
// each kernel block needs because several glibc files reuse constant names
// (hp0/hp1/big/mp1/mp2/split/T/N/A/OFF) with different values, and (4) the
// FMA/TOINT/errno/FP-exception branches the qualified build does not
// compile are omitted or #if'd out exactly where glibc's own guards sit.

#include "hostlineage_internal.hpp"

#include <kalos/detail/hostlineage.hpp>

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace kalos::detail::hostlineage_src {

// glibc targets this transcription was proven against evaluate double_t as
// plain double (FLT_EVAL_METHOD == 0); the assert below keeps that identity
// load-bearing rather than assumed.
static_assert(FLT_EVAL_METHOD == 0,
              "hostlineage kernels require a target with FLT_EVAL_METHOD == 0");

namespace {

// ---------------------------------------------------------------------------
// Plumbing substitutions for glibc's math_private.h / math_config.h /
// math-barriers.h / fenv_private.h / math-underflow.h (bit-equivalent only;
// see public header for the semantics decisions).
// ---------------------------------------------------------------------------

// glibc's asuint64/asdouble use union type-punning; memcpy is the same
// operation, UB-clean, and optimized identically by GCC/Clang at -O1+.
inline std::uint64_t asuint64(double f) {
  std::uint64_t i;
  std::memcpy(&i, &f, sizeof i);
  return i;
}
inline double asdouble(std::uint64_t i) {
  double f;
  std::memcpy(&f, &i, sizeof f);
  return f;
}

// math_config.h wants_rounding configuration: return-value bits only.
constexpr int WANT_ROUNDING = 1;

// math_narrow_eval is the identity under FLT_EVAL_METHOD == 0 (asserted).
inline double math_narrow_eval(double x) { return x; }

// math-barriers.h: only needed to defeat speculative code motion around FP
// exception signalling, which this transcription does not perform.
template <typename T>
inline T math_opt_barrier(T x) { return x; }

// Branch hints (cdefs.h). glibc's own cdefs.h may already define these as
// __builtin_expect wrappers (pure hints, no effect on bits); undef first.
#ifdef __glibc_unlikely
#undef __glibc_unlikely
#endif
#ifdef __glibc_likely
#undef __glibc_likely
#endif
#define __glibc_unlikely(cond) (cond)
#define __glibc_likely(cond) (cond)

// s_math_uflow.c / s_math_oflow.c return-value bits (nearest rounding).
inline double __math_uflow(std::uint32_t sign) { return sign ? -0.0 : 0.0; }
inline double __math_oflow(std::uint32_t sign) {
  return sign ? -INFINITY : INFINITY;
}
// s_math_divzero.c: ±1/0.
inline double __math_divzero(std::uint32_t sign) {
  return sign ? -INFINITY : INFINITY;
}
// s_math_invalid.c: default-quieting NaN.
inline double __math_invalid(double x) { return (x - x) / (x - x); }

// math_config.h check_oflow/check_uflow with WANT_ERRNO = 0: pass-throughs.
inline double check_oflow(double x) { return x; }
inline double check_uflow(double x) { return x; }

// math_config.h issignaling_inline, LITTLE branch (HIGH_ORDER_BIT_IS_SET_FOR_
// SNAN = 0 on x86-64).
inline int issignaling_inline(double x) {
  std::uint64_t ix = asuint64(x);
  return 2 * (ix ^ 0x0008000000000000) > 2 * 0x7ff8000000000000ULL;
}

// ---------------------------------------------------------------------------
// dla.h double-length arithmetic (non-FMA branch: DLA_FMS is never defined
// on the qualified target — __FP_FAST_FMA absent). glibc implements these as
// statement-expression macros; they are verbatim reference-output functions
// here so the call sites below are character-identical to glibc.
// ---------------------------------------------------------------------------

#define CN 134217729.0 /* 1+2**27, IEEE double 0x41a0000002000000 */

// Exact addition of two single-length floating point numbers, Dekker.
inline void EADD(double x, double y, double& z, double& zz) {
  z = (x) + (y);
  zz = (std::fabs(x) > std::fabs(y)) ? (((x) - (z)) + (y)) : (((y) - (z)) + (x));
}

// Exact subtraction of two single-length floating point numbers, Dekker.
inline void ESUB(double x, double y, double& z, double& zz) {
  z = (x) - (y);
  zz = (std::fabs(x) > std::fabs(y)) ? (((x) - (z)) - (y)) : ((x) - ((y) + (z)));
}

// Exact multiplication of two single-length floating point numbers, Veltkamp.
inline void EMULV(double x, double y, double& z, double& zz) {
  double p, hx, tx, hy, ty;
  p = CN * (x);
  hx = ((x) - p) + p;
  tx = (x) - hx;
  p = CN * (y);
  hy = ((y) - p) + p;
  ty = (y) - hy;
  z = (x) * (y);
  zz = (((hx * hy - z) + hx * ty) + tx * hy) + tx * ty;
}

// ===========================================================================
// sin / cos — s_sin.c + branred.c (IBM APML)
// ===========================================================================
namespace sincos {

// usncs.h constants (LITTLE_ENDI-independent block). Unused members of
// usncs.h (aa, bb, mp3) are omitted: no vendored code path reads them.
static const double s1 = -0x1.5555555555555p-3;   /* -0.16666666666666666     */
static const double s2 = 0x1.1111111110ECEp-7;    /*  0.0083333333333323288   */
static const double s3 = -0x1.A01A019DB08B8p-13;  /* -0.00019841269834414642  */
static const double s4 = 0x1.71DE27B9A7ED9p-19;   /*  2.755729806860771e-06   */
static const double s5 = -0x1.ADDFFC2FCDF59p-26;  /* -2.5022014848318398e-08  */
static const double big = 0x1.8000000000000p45;   /*  52776558133248          */
static const double hp0 = 0x1.921FB54442D18p0;    /*  1.5707963267948966      */
static const double hp1 = 0x1.1A62633145C07p-54;  /*  6.123233995736766e-17   */
static const double mp1 = 0x1.921FB58000000p0;    /*  1.5707963407039642      */
static const double mp2 = -0x1.DDE973C000000p-27; /* -1.3909067564377153e-08  */
static const double pp3 = -0x1.CB3B398000000p-55; /* -4.9789962314799099e-17  */
static const double pp4 = -0x1.d747f23e32ed7p-83; /* -1.9034889620193266e-25  */
static const double hpinv = 0x1.45F306DC9C883p-1; /*  0.63661977236758138     */
static const double toint = 0x1.8000000000000p52; /*  6755399441055744        */

// branred.h constants, LITTLE_ENDI branch verbatim (branred's mp2 differs in
// its low word from usncs.h mp2 — both are kept exactly as glibc ships them).
namespace branred_c {

static const mynumber
/**/           t576 = {{0x00000000, 0x63f00000}},  /* 2 ^ 576  */
/**/          tm600 = {{0x00000000, 0x1a700000}},  /* 2 ^- 600 */
/**/           tm24 = {{0x00000000, 0x3e700000}},  /* 2 ^- 24  */
/**/            big = {{0x00000000, 0x43380000}},  /*  6755399441055744      */
/**/           big1 = {{0x00000000, 0x43580000}},  /* 27021597764222976      */
/**/            hp0 = {{0x54442D18, 0x3FF921FB}},  /* 1.5707963267948966     */
/**/            hp1 = {{0x33145C07, 0x3C91A626}},  /* 6.123233995736766e-17  */
/**/            mp1 = {{0x58000000, 0x3FF921FB}},  /* 1.5707963407039642     */
/**/            mp2 = {{0x40000000, IW(0xBE4DDE97)}};  /*-1.3909067675399456e-08 */

static const double toverp[75] = { /*  2/ PI base 24*/
  10680707.0,  7228996.0,  1387004.0,  2578385.0, 16069853.0,
  12639074.0,  9804092.0,  4427841.0, 16666979.0, 11263675.0,
  12935607.0,  2387514.0,  4345298.0, 14681673.0,  3074569.0,
  13734428.0, 16653803.0,  1880361.0, 10960616.0,  8533493.0,
   3062596.0,  8710556.0,  7349940.0,  6258241.0,  3772886.0,
   3769171.0,  3798172.0,  8675211.0, 12450088.0,  3874808.0,
   9961438.0,   366607.0, 15675153.0,  9132554.0,  7151469.0,
   3571407.0,  2607881.0, 12013382.0,  4155038.0,  6285869.0,
   7677882.0, 13102053.0, 15825725.0,   473591.0,  9065106.0,
  15363067.0,  6271263.0,  9264392.0,  5636912.0,  4652155.0,
   7056368.0, 13614112.0, 10155062.0,  1944035.0,  9527646.0,
  15080200.0,  6658437.0,  6231200.0,  6832269.0, 16767104.0,
   5075751.0,  3212806.0,  1398474.0,  7579849.0,  6349435.0,
  12618859.0,  4703257.0, 12806093.0, 14477321.0,  2786137.0,
  12875403.0,  9837734.0, 14528324.0, 13719321.0,   343717.0 };

static const double split = CN;	/* 2^27 + 1 */

// branred.c __branred. The math_private.h include and SECTION plumbing are
// the only omissions; the body is character-faithful (int4 -> int, mynumber
// union punning identical to glibc).
int branred(double x, double* a, double* aa) {
  int i, k;
  mynumber u, gor;
  double r[6], s, t, sum, b, bb, sum1, sum2, b1, bb1, b2, bb2, x1, x2, t1, t2;

  x *= tm600.x;
  t = x * split; /* split x to two numbers */
  x1 = t - (t - x);
  x2 = x - x1;
  sum = 0;
  u.x = x1;
  k = (u.i[HIGH_HALF] >> 20) & 2047;
  k = (k - 450) / 24;
  if (k < 0)
    k = 0;
  gor.x = t576.x;
  gor.i[HIGH_HALF] -= ((k * 24) << 20);
  for (i = 0; i < 6; i++) {
    r[i] = x1 * toverp[k + i] * gor.x;
    gor.x *= tm24.x;
  }
  for (i = 0; i < 3; i++) {
    s = (r[i] + big.x) - big.x;
    sum += s;
    r[i] -= s;
  }
  t = 0;
  for (i = 0; i < 6; i++)
    t += r[5 - i];
  bb = (((((r[0] - t) + r[1]) + r[2]) + r[3]) + r[4]) + r[5];
  s = (t + big.x) - big.x;
  sum += s;
  t -= s;
  b = t + bb;
  bb = (t - b) + bb;
  s = (sum + big1.x) - big1.x;
  sum -= s;
  b1 = b;
  bb1 = bb;
  sum1 = sum;
  sum = 0;

  u.x = x2;
  k = (u.i[HIGH_HALF] >> 20) & 2047;
  k = (k - 450) / 24;
  if (k < 0)
    k = 0;
  gor.x = t576.x;
  gor.i[HIGH_HALF] -= ((k * 24) << 20);
  for (i = 0; i < 6; i++) {
    r[i] = x2 * toverp[k + i] * gor.x;
    gor.x *= tm24.x;
  }
  for (i = 0; i < 3; i++) {
    s = (r[i] + big.x) - big.x;
    sum += s;
    r[i] -= s;
  }
  t = 0;
  for (i = 0; i < 6; i++)
    t += r[5 - i];
  bb = (((((r[0] - t) + r[1]) + r[2]) + r[3]) + r[4]) + r[5];
  s = (t + big.x) - big.x;
  sum += s;
  t -= s;
  b = t + bb;
  bb = (t - b) + bb;
  s = (sum + big1.x) - big1.x;
  sum -= s;

  b2 = b;
  bb2 = bb;
  sum2 = sum;

  sum = sum1 + sum2;
  b = b1 + b2;
  bb = (std::fabs(b1) > std::fabs(b2)) ? (b1 - b) + b2 : (b2 - b) + b1;
  if (b > 0.5) {
    b -= 1.0;
    sum += 1.0;
  } else if (b < -0.5) {
    b += 1.0;
    sum -= 1.0;
  }
  s = b + (bb + bb1 + bb2);
  t = ((b - s) + bb) + (bb1 + bb2);
  b = s * split;
  t1 = b - (b - s);
  t2 = s - t1;
  b = s * hp0.x;
  bb = ((((t1 * mp1.x - b) + t1 * mp2.x) + t2 * mp1.x) +
        (t2 * mp2.x + s * hp1.x + t * hp0.x));
  s = b + bb;
  t = (b - s) + bb;
  *a = s;
  *aa = t;
  return ((int)sum) & 3; /* return quarter of unit circle */
}

}  // namespace branred_c

// s_sin.c helper macros -> verbatim inline functions (-Wpedantic: statement
// expressions are not used; call sites keep glibc's syntax).
#define POLYNOMIAL2(xx) ((((s5 * (xx) + s4) * (xx) + s3) * (xx) + s2) * (xx))
#define POLYNOMIAL(xx) (POLYNOMIAL2(xx) + s1)

inline double TAYLOR_SIN(double xx, double x, double dx) {
  double t = ((POLYNOMIAL(xx) * (x) - 0.5 * (dx)) * (xx) + (dx));
  double res = (x) + t;
  return res;
}

inline void SINCOS_TABLE_LOOKUP(const mynumber& u, double& sn, double& ssn,
                                double& cs, double& ccs) {
  int4 k = u.i[LOW_HALF] << 2;
  sn = sincostab.x[k];
  ssn = sincostab.x[k + 1];
  cs = sincostab.x[k + 2];
  ccs = sincostab.x[k + 3];
}

static const double
  sn3 = -1.66666666666664880952546298448555E-01,
  sn5 = 8.33333214285722277379541354343671E-03,
  cs2 = 4.99999999999999999999950396842453E-01,
  cs4 = -4.16666666666664434524222570944589E-02,
  cs6 = 1.38888874007937613028114285595617E-03;

inline double do_cos(double x, double dx) {
  mynumber u;

  if (x < 0)
    dx = -dx;

  u.x = big + std::fabs(x);
  x = std::fabs(x) - (u.x - big) + dx;

  double xx, s, sn, ssn, c, cs, ccs, cor;
  xx = x * x;
  s = x + x * xx * (sn3 + xx * sn5);
  c = xx * (cs2 + xx * (cs4 + xx * cs6));
  SINCOS_TABLE_LOOKUP(u, sn, ssn, cs, ccs);
  cor = (ccs - s * ssn - cs * c) - sn * s;
  return cs + cor;
}

inline double do_sin(double x, double dx) {
  double xold = x;
  /* Max ULP is 0.501 if |x| < 0.126, otherwise ULP is 0.518.  */
  if (std::fabs(x) < 0.126)
    return TAYLOR_SIN(x * x, x, dx);

  mynumber u;

  if (x <= 0)
    dx = -dx;
  u.x = big + std::fabs(x);
  x = std::fabs(x) - (u.x - big);

  double xx, s, sn, ssn, c, cs, ccs, cor;
  xx = x * x;
  s = x + (dx + x * xx * (sn3 + xx * sn5));
  c = x * dx + xx * (cs2 + xx * (cs4 + xx * cs6));
  SINCOS_TABLE_LOOKUP(u, sn, ssn, cs, ccs);
  cor = (ssn + s * ccs - sn * c) + cs * s;
  return std::copysign(sn + cor, xold);
}

inline int4 reduce_sincos(double x, double* a, double* da) {
  mynumber v;

  double t = (x * hpinv + toint);
  double xn = t - toint;
  v.x = t;
  double y = (x - xn * mp1) - xn * mp2;
  int4 n = v.i[LOW_HALF] & 3;

  double b, db, t1, t2;
  t1 = xn * pp3;
  t2 = y - t1;
  db = (y - t2) - t1;

  t1 = xn * pp4;
  b = t2 - t1;
  db += (t2 - b) - t1;

  *a = b;
  *da = db;
  return n;
}

inline double do_sincos(double a, double da, int4 n) {
  double retval;

  if (n & 1)
    retval = do_cos(a, da);
  else
    retval = do_sin(a, da);

  return (n & 2) ? -retval : retval;
}

// __sin / __cos. SET_RESTORE_ROUND_53BIT(FE_TONEAREST) is a no-op here
// (round-to-nearest is the one mode the entire numeric surface assumes), and
// math_check_force_underflow/__set_errno cannot affect returned bits.
double sin_kernel(double x) {
  double t, a, da;
  mynumber u;
  int4 k, m, n;
  double retval = 0;

  u.x = x;
  m = u.i[HIGH_HALF];
  k = 0x7fffffff & m; /* no sign           */
  if (k < 0x3e500000) /* if x->0 =>sin(x)=x */
  {
    retval = x;
  }
  /*--------------------------- 2^-26<|x|< 0.855469---------------------- */
  else if (k < 0x3feb6000) {
    retval = do_sin(x, 0);
  }

  /*----------------------- 0.855469  <|x|<2.426265  ----------------------*/
  else if (k < 0x400368fd) {
    t = hp0 - std::fabs(x);
    retval = std::copysign(do_cos(t, hp1), x);
  }

  /*-------------------------- 2.426265<|x|< 105414350 ----------------------*/
  else if (k < 0x419921FB) {
    n = reduce_sincos(x, &a, &da);
    retval = do_sincos(a, da, n);
  }

  /* --------------------105414350 <|x| <2^1024------------------------------*/
  else if (k < 0x7ff00000) {
    n = branred_c::branred(x, &a, &da);
    retval = do_sincos(a, da, n);
  }
  /*--------------------- |x| > 2^1024 ----------------------------------*/
  else {
    retval = x / x;
  }

  return retval;
}

double cos_kernel(double x) {
  double y, a, da;
  mynumber u;
  int4 k, m, n;

  double retval = 0;

  u.x = x;
  m = u.i[HIGH_HALF];
  k = 0x7fffffff & m;

  /* |x|<2^-27 => cos(x)=1 */
  if (k < 0x3e400000)
    retval = 1.0;

  else if (k < 0x3feb6000) { /* 2^-27 < |x| < 0.855469 */
    retval = do_cos(x, 0);
  }

  else if (k < 0x400368fd) { /* 0.855469  <|x|<2.426265  */
    y = hp0 - std::fabs(x);
    a = y + hp1;
    da = (y - a) + hp1;
    retval = do_sin(a, da);
  }

  else if (k < 0x419921FB) { /* 2.426265<|x|< 105414350 */
    n = reduce_sincos(x, &a, &da);
    retval = do_sincos(a, da, n + 1);
  }

  /* 105414350 <|x| <2^1024 */
  else if (k < 0x7ff00000) {
    n = branred_c::branred(x, &a, &da);
    retval = do_sincos(a, da, n + 1);
  }

  else {
    retval = x / x; /* |x| > 2^1024 */
  }

  return retval;
}

}  // namespace sincos

// ===========================================================================
// exp — e_exp.c (ARM table-driven). TOINT_INTRINSICS = 0 (Shift-path; glibc
// never sets TOINT_INTRINSICS on x86-64) and the non-FMA code are the only
// compiled variants, matching the qualified libm (public header: no-IFUNC
// audit).
// ===========================================================================
namespace exp {

constexpr int N = 1 << 7;  // EXP_TABLE_BITS
// C2..C5 == poly[5-EXP_POLY_ORDER .. 8-EXP_POLY_ORDER] == poly[0..3].

inline std::uint32_t top12(double x) { return asuint64(x) >> 52; }

// specialcase (e_exp.c, WANT_ERRNO=0; the math_force_eval underflow raise
// cannot change returned bits and is omitted).
inline double specialcase(double tmp, std::uint64_t sbits, std::uint64_t ki) {
  double scale, y;

  if ((ki & 0x80000000) == 0) {
    /* k > 0, the exponent of scale might have overflowed by <= 460.  */
    sbits -= 1009ull << 52;
    scale = asdouble(sbits);
    y = 0x1p1009 * (scale + scale * tmp);
    return check_oflow(y);
  }
  /* k < 0, need special care in the subnormal range.  */
  sbits += 1022ull << 52;
  scale = asdouble(sbits);
  y = scale + scale * tmp;
  if (y < 1.0) {
    double hi, lo;
    lo = scale - y + scale * tmp;
    hi = 1.0 + y;
    lo = 1.0 - hi + y + lo;
    y = math_narrow_eval(hi + lo) - 1.0;
    /* Avoid -0.0 with downward rounding.  */
    if (WANT_ROUNDING && y == 0.0)
      y = 0.0;
  }
  y = 0x1p-1022 * y;
  return check_uflow(y);
}

double exp_kernel(double x) {
  std::uint32_t abstop;
  std::uint64_t ki, idx, top, sbits;
  double kd, z, r, r2, scale, tail, tmp;

  abstop = top12(x) & 0x7ff;
  if (__glibc_unlikely(abstop - top12(0x1p-54) >=
                       top12(512.0) - top12(0x1p-54))) {
    if (abstop - top12(0x1p-54) >= 0x80000000)
      /* Avoid spurious underflow for tiny x.  */
      /* Note: 0 is common input.  */
      return WANT_ROUNDING ? 1.0 + x : 1.0;
    if (abstop >= top12(1024.0)) {
      if (asuint64(x) == asuint64(-INFINITY))
        return 0.0;
      if (abstop >= top12(INFINITY))
        return 1.0 + x;
      if (asuint64(x) >> 63)
        return __math_uflow(0);
      else
        return __math_oflow(0);
    }
    /* Large x is special cased below.  */
    abstop = 0;
  }

  /* exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].  */
  /* x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].  */
  z = exp_data.invln2N * x;
  /* TOINT_INTRINSICS == 0 branch (glibc x86-64).  */
  kd = math_narrow_eval(z + exp_data.shift);
  ki = asuint64(kd);
  kd -= exp_data.shift;
  r = x + kd * exp_data.negln2hiN + kd * exp_data.negln2loN;
  /* 2^(k/N) ~= scale * (1 + tail).  */
  idx = 2 * (ki % N);
  top = ki << (52 - 7);  // EXP_TABLE_BITS
  tail = asdouble(exp_data.tab[idx]);
  /* This is only a valid scale when -1023*N < k < 1024*N.  */
  sbits = exp_data.tab[idx + 1] + top;
  /* exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).  */
  r2 = r * r;
  tmp = tail + r + r2 * (exp_data.poly[0] + r * exp_data.poly[1]) +
        r2 * r2 * (exp_data.poly[2] + r * exp_data.poly[3]);
  if (__glibc_unlikely(abstop == 0))
    return specialcase(tmp, sbits, ki);
  scale = asdouble(sbits);
  return scale + scale * tmp;
}

}  // namespace exp

// ===========================================================================
// log — e_log.c (ARM table-driven, non-FMA branch: r built from tab2).
// ===========================================================================
namespace log {

constexpr std::uint64_t OFF = 0x3fe6000000000000ULL;
constexpr int N = 1 << 7;  // LOG_TABLE_BITS

// top16: top 16 bits of a double.
inline std::uint32_t top16(double x) { return asuint64(x) >> 48; }

double log_kernel(double x) {
  double w, z, r, r2, r3, y, invc, logc, kd, hi, lo;
  std::uint64_t ix, iz, tmp;
  std::uint32_t top;
  int k, i;

  ix = asuint64(x);
  top = top16(x);

  const std::uint64_t LO = asuint64(1.0 - 0x1p-4);
  const std::uint64_t HI = asuint64(1.0 + 0x1.09p-4);
  if (__glibc_unlikely(ix - LO < HI - LO)) {
    /* Handle close to 1.0 inputs separately.  */
    /* Fix sign of zero with downward rounding when x==1.  */
    if (WANT_ROUNDING && __glibc_unlikely(ix == asuint64(1.0)))
      return 0;
    r = x - 1.0;
    r2 = r * r;
    r3 = r * r2;
    y = r3 * (log_data.poly1[1] + r * log_data.poly1[2] +
               r2 * log_data.poly1[3] +
               r3 * (log_data.poly1[4] + r * log_data.poly1[5] +
                      r2 * log_data.poly1[6] +
                      r3 * (log_data.poly1[7] + r * log_data.poly1[8] +
                             r2 * log_data.poly1[9] + r3 * log_data.poly1[10])));
    /* Worst-case error is around 0.507 ULP.  */
    w = r * 0x1p27;
    double rhi = r + w - w;
    double rlo = r - rhi;
    w = rhi * rhi * log_data.poly1[0]; /* B[0] == -0.5.  */
    hi = r + w;
    lo = r - hi + w;
    lo += log_data.poly1[0] * rlo * (rhi + r);
    y += lo;
    y += hi;
    return y;
  }
  if (__glibc_unlikely(top - 0x0010 >= 0x7ff0 - 0x0010)) {
    /* x < 0x1p-1022 or inf or nan.  */
    if (ix * 2 == 0)
      return __math_divzero(1);
    if (ix == asuint64(INFINITY)) /* log(inf) == inf.  */
      return x;
    if ((top & 0x8000) || (top & 0x7ff0) == 0x7ff0)
      return __math_invalid(x);
    /* x is subnormal, normalize it.  */
    ix = asuint64(x * 0x1p52);
    ix -= 52ULL << 52;
  }

  /* x = 2^k z; where z is in range [OFF,2*OFF) and exact.  */
  tmp = ix - OFF;
  i = (tmp >> (52 - 7)) % N;  // LOG_TABLE_BITS
  k = (std::int64_t)tmp >> 52; /* arithmetic shift */
  iz = ix - (tmp & 0xfffULL << 52);
  invc = log_data.tab[i].invc;
  logc = log_data.tab[i].logc;
  z = asdouble(iz);

  /* log(x) = log1p(z/c-1) + log(c) + k*Ln2.  */
  /* non-FMA branch (rounding error: 0x1p-55/N + 0x1p-66).  */
  r = (z - log_data.tab2[i].chi - log_data.tab2[i].clo) * invc;
  kd = (double)k;

  /* hi + lo = r + log(c) + k*Ln2.  */
  w = kd * log_data.ln2hi + logc;
  hi = w + r;
  lo = w - hi + r + kd * log_data.ln2lo;

  /* log(x) = lo + (log1p(r) - r) + hi.  */
  r2 = r * r; /* rounding error: 0x1p-54/N^2.  */
  y = lo + r2 * log_data.poly[0] +
      r * r2 * (log_data.poly[1] + r * log_data.poly[2] +
                r2 * (log_data.poly[3] + r * log_data.poly[4])) +
      hi;
  return y;
}

}  // namespace log

// ===========================================================================
// pow — e_pow.c (ARM table-driven; own log_inline/exp_inline variants,
// non-FMA branches, Shift-path reduction, WANT_ERRNO=0, WANT_ROUNDING=1).
// ===========================================================================
namespace pow {

constexpr std::uint64_t OFF = 0x3fe6955500000000ULL;
constexpr int LOG_N = 1 << 7;      // POW_LOG_TABLE_BITS
constexpr int EXP_N = 1 << 7;      // EXP_TABLE_BITS
constexpr std::uint32_t SIGN_BIAS = 0x800 << 7;

inline std::uint32_t top12(double x) { return asuint64(x) >> 52; }

// log_inline: y + *tail = log(x); IX normalized in the subnormal range with
// the sign bit used for the exponent.
inline double log_inline(std::uint64_t ix, double* tail) {
  double z, r, y, invc, logc, logctail, kd, hi, t1, t2, lo, lo1, lo2, p;
  std::uint64_t iz, tmp;
  int k, i;

  /* x = 2^k z; where z is in range [OFF,2*OFF) and exact.  */
  tmp = ix - OFF;
  i = (tmp >> (52 - 7)) % LOG_N;  // POW_LOG_TABLE_BITS
  k = (std::int64_t)tmp >> 52; /* arithmetic shift */
  iz = ix - (tmp & 0xfffULL << 52);
  z = asdouble(iz);
  kd = (double)k;

  /* log(x) = k*Ln2 + log(c) + log1p(z/c-1).  */
  invc = pow_log_data.tab[i].invc;
  logc = pow_log_data.tab[i].logc;
  logctail = pow_log_data.tab[i].logctail;

  /* non-FMA branch: split z so rhi, rlo and rhi*rhi are exact.  */
  double zhi = asdouble((iz + (1ULL << 31)) & (-1ULL << 32));
  double zlo = z - zhi;
  double rhi = zhi * invc - 1.0;
  double rlo = zlo * invc;
  r = rhi + rlo;

  /* k*Ln2 + log(c) + r.  */
  t1 = kd * pow_log_data.ln2hi + logc;
  t2 = t1 + r;
  lo1 = kd * pow_log_data.ln2lo + logctail;
  lo2 = t1 - t2 + r;

  double ar, ar2, ar3, lo3, lo4;
  ar = pow_log_data.poly[0] * r; /* A[0] = -0.5.  */
  ar2 = r * ar;
  ar3 = r * ar2;
  /* k*Ln2 + log(c) + r + A[0]*r*r.  */
  /* non-FMA branch.  */
  double arhi = pow_log_data.poly[0] * rhi;
  double arhi2 = rhi * arhi;
  hi = t2 + arhi2;
  lo3 = rlo * (ar + arhi);
  lo4 = t2 - hi + arhi2;
  /* p = log1p(r) - r - A[0]*r*r.  */
  p = (ar3 * (pow_log_data.poly[1] + r * pow_log_data.poly[2] +
              ar2 * (pow_log_data.poly[3] + r * pow_log_data.poly[4] +
                     ar2 * (pow_log_data.poly[5] + r * pow_log_data.poly[6]))));
  lo = lo1 + lo2 + lo3 + lo4 + p;
  y = hi + lo;
  *tail = hi - y + lo;
  return y;
}

// specialcase (e_pow.c; signed-scale variant; WANT_ERRNO=0).
inline double specialcase(double tmp, std::uint64_t sbits, std::uint64_t ki) {
  double scale, y;

  if ((ki & 0x80000000) == 0) {
    /* k > 0, the exponent of scale might have overflowed by <= 460.  */
    sbits -= 1009ull << 52;
    scale = asdouble(sbits);
    y = 0x1p1009 * (scale + scale * tmp);
    return check_oflow(y);
  }
  /* k < 0, need special care in the subnormal range.  */
  sbits += 1022ull << 52;
  /* Note: sbits is signed scale.  */
  scale = asdouble(sbits);
  y = scale + scale * tmp;
  if (std::fabs(y) < 1.0) {
    double hi, lo, one = 1.0;
    if (y < 0.0)
      one = -1.0;
    lo = scale - y + scale * tmp;
    hi = one + y;
    lo = one - hi + y + lo;
    y = math_narrow_eval(hi + lo) - one;
    /* Fix the sign of 0.  */
    if (y == 0.0)
      y = asdouble(sbits & 0x8000000000000000);
  }
  y = 0x1p-1022 * y;
  return check_uflow(y);
}

// exp_inline: sign_bias*exp(x + xtail), Shift-path reduction (non-FMA).
inline double exp_inline(double x, double xtail, std::uint32_t sign_bias) {
  std::uint32_t abstop;
  std::uint64_t ki, idx, top, sbits;
  double kd, z, r, r2, scale, tail, tmp;

  abstop = top12(x) & 0x7ff;
  if (__glibc_unlikely(abstop - top12(0x1p-54) >=
                       top12(512.0) - top12(0x1p-54))) {
    if (abstop - top12(0x1p-54) >= 0x80000000) {
      /* Avoid spurious underflow for tiny x.  */
      /* Note: 0 is common input.  */
      double one = WANT_ROUNDING ? 1.0 + x : 1.0;
      return sign_bias ? -one : one;
    }
    if (abstop >= top12(1024.0)) {
      /* Note: inf and nan are already handled.  */
      if (asuint64(x) >> 63)
        return __math_uflow(sign_bias);
      else
        return __math_oflow(sign_bias);
    }
    /* Large x is special cased below.  */
    abstop = 0;
  }

  z = exp_data.invln2N * x;
  /* TOINT_INTRINSICS == 0 branch.  */
  kd = math_narrow_eval(z + exp_data.shift);
  ki = asuint64(kd);
  kd -= exp_data.shift;
  r = x + kd * exp_data.negln2hiN + kd * exp_data.negln2loN;
  /* The code assumes 2^-200 < |xtail| < 2^-8/N.  */
  r += xtail;
  idx = 2 * (ki % EXP_N);
  top = (ki + sign_bias) << (52 - 7);  // EXP_TABLE_BITS
  tail = asdouble(exp_data.tab[idx]);
  sbits = exp_data.tab[idx + 1] + top;
  r2 = r * r;
  tmp = tail + r + r2 * (exp_data.poly[0] + r * exp_data.poly[1]) +
        r2 * r2 * (exp_data.poly[2] + r * exp_data.poly[3]);
  if (__glibc_unlikely(abstop == 0))
    return specialcase(tmp, sbits, ki);
  scale = asdouble(sbits);
  return scale + scale * tmp;
}

// Returns 0 if not int, 1 if odd int, 2 if even int.
inline int checkint(std::uint64_t iy) {
  int e = iy >> 52 & 0x7ff;
  if (e < 0x3ff)
    return 0;
  if (e > 0x3ff + 52)
    return 2;
  if (iy & ((1ULL << (0x3ff + 52 - e)) - 1))
    return 0;
  if (iy & (1ULL << (0x3ff + 52 - e)))
    return 1;
  return 2;
}

// Returns 1 if input is the bit representation of 0, infinity or nan.
inline int zeroinfnan(std::uint64_t i) {
  return 2 * i - 1 >= 2 * asuint64(INFINITY) - 1;
}

double pow_kernel(double x, double y) {
  std::uint32_t sign_bias = 0;
  std::uint64_t ix, iy;
  std::uint32_t topx, topy;

  ix = asuint64(x);
  iy = asuint64(y);
  topx = top12(x);
  topy = top12(y);
  if (__glibc_unlikely(topx - 0x001 >= 0x7ff - 0x001 ||
                       (topy & 0x7ff) - 0x3be >= 0x43e - 0x3be)) {
    /* Special cases: (x < 0x1p-126 or inf or nan) or
       (|y| < 0x1p-65 or |y| >= 0x1p63 or nan).  */
    if (__glibc_unlikely(zeroinfnan(iy))) {
      if (2 * iy == 0)
        return issignaling_inline(x) ? x + y : 1.0;
      if (ix == asuint64(1.0))
        return issignaling_inline(y) ? x + y : 1.0;
      if (2 * ix > 2 * asuint64(INFINITY) || 2 * iy > 2 * asuint64(INFINITY))
        return x + y;
      if (2 * ix == 2 * asuint64(1.0))
        return 1.0;
      if ((2 * ix < 2 * asuint64(1.0)) == !(iy >> 63))
        return 0.0; /* |x|<1 && y==inf or |x|>1 && y==-inf.  */
      return y * y;
    }
    if (__glibc_unlikely(zeroinfnan(ix))) {
      double x2 = x * x;
      if (ix >> 63 && checkint(iy) == 1) {
        x2 = -x2;
        sign_bias = 1;
      }
      /* WANT_ERRNO == 1 (glibc default — this branch changes RETURN BITS,
         not just errno: pow(-0, odd-negative) must be -inf, while the plain
         1/x2 path below yields +0-divide): kept verbatim.  */
      if (2 * ix == 0 && iy >> 63)
        return __math_divzero(sign_bias);
      return iy >> 63 ? math_opt_barrier(1 / x2) : x2;
    }
    /* Here x and y are non-zero finite.  */
    if (ix >> 63) {
      /* Finite x < 0.  */
      int yint = checkint(iy);
      if (yint == 0)
        return __math_invalid(x);
      if (yint == 1)
        sign_bias = SIGN_BIAS;
      ix &= 0x7fffffffffffffff;
      topx &= 0x7ff;
    }
    if ((topy & 0x7ff) - 0x3be >= 0x43e - 0x3be) {
      /* Note: sign_bias == 0 here because y is not odd.  */
      if (ix == asuint64(1.0))
        return 1.0;
      if ((topy & 0x7ff) < 0x3be) {
        /* |y| < 2^-65, x^y ~= 1 + y*log(x).  */
        if (WANT_ROUNDING)
          return ix > asuint64(1.0) ? 1.0 + y : 1.0 - y;
        else
          return 1.0;
      }
      return (ix > asuint64(1.0)) == (topy < 0x800) ? __math_oflow(0)
                                                    : __math_uflow(0);
    }
    if (topx == 0) {
      /* Normalize subnormal x so exponent becomes negative.  */
      ix = asuint64(x * 0x1p52);
      ix &= 0x7fffffffffffffff;
      ix -= 52ULL << 52;
    }
  }

  double lo;
  double hi = log_inline(ix, &lo);
  double ehi, elo;
  /* non-FMA branch.  */
  double yhi = asdouble(iy & -1ULL << 27);
  double ylo = y - yhi;
  double lhi = asdouble(asuint64(hi) & -1ULL << 27);
  double llo = hi - lhi + lo;
  ehi = yhi * lhi;
  elo = ylo * lhi + y * llo; /* |elo| < |ehi| * 2^-25.  */
  return exp_inline(ehi, elo, sign_bias);
}

}  // namespace pow

// ===========================================================================
// hypot — e_hypot.c (Borges MyHypot3, non-FMA correction branch; errno and
// the underflow-exception raise cannot change returned bits).
// ===========================================================================
namespace hypotk {

constexpr double SCALE = 0x1p-600;
constexpr double LARGE_VAL = 0x1p+511;
constexpr double TINY_VAL = 0x1p-459;
constexpr double EPS = 0x1p-54;

// kernel: inputs adjusted so ax >= ay >= 0 and squaring ax, ay and
// (ax - ay) neither overflows nor underflows.
inline double kernel(double ax, double ay) {
  double t1, t2;
  /* non-FMA branch (Borges correction).  */
  double h = std::sqrt(ax * ax + ay * ay);
  if (h <= 2.0 * ay) {
    double delta = h - ay;
    t1 = ax * (2.0 * delta - ax);
    t2 = (delta - 2.0 * (ax - ay)) * delta;
  } else {
    double delta = h - ax;
    t1 = 2.0 * delta * (ax - 2.0 * ay);
    t2 = (4.0 * delta - ay) * ay + delta * delta;
  }

  h -= (t1 + t2) / (2.0 * h);
  return h;
}

double hypot_kernel(double x, double y) {
  if (!std::isfinite(x) || !std::isfinite(y)) {
    if ((std::isinf(x) || std::isinf(y)) && !issignaling_inline(x) &&
        !issignaling_inline(y))
      return INFINITY;
    return x + y;
  }

  x = std::fabs(x);
  y = std::fabs(y);

  // USE_FMAX_BUILTIN branch of glibc is bit-identical to the ternary for the
  // finite, non-negative, NaN-free inputs that survive the prologue above.
  double ax = x < y ? y : x;
  double ay = x < y ? x : y;

  /* If ax is huge, scale both inputs down.  */
  if (__glibc_unlikely(ax > LARGE_VAL)) {
    if (__glibc_unlikely(ay <= ax * EPS))
      return math_narrow_eval(ax + ay);

    return math_narrow_eval(kernel(ax * SCALE, ay * SCALE) / SCALE);
  }

  /* If ay is tiny, scale both inputs up.  */
  if (__glibc_unlikely(ay < TINY_VAL)) {
    if (__glibc_unlikely(ax >= ay / EPS))
      return math_narrow_eval(ax + ay);

    ax = math_narrow_eval(kernel(ax / SCALE, ay / SCALE) * SCALE);
    return ax;
  }

  /* Common case: ax is not huge and ay is not tiny.  */
  if (__glibc_unlikely(ay <= ax * EPS))
    return ax + ay;

  return kernel(ax, ay);
}

}  // namespace hypotk

// ===========================================================================
// atan2 — e_atan2.c (IBM APML uatan2 + dla.h Dekker arithmetic).
// SET_RESTORE_ROUND(FE_TONEAREST) is a no-op; the DBL_MIN force-underflow
// raise cannot change returned bits.
// ===========================================================================
namespace atan2 {

// atnat2.h constants, LITTLE_ENDI branch (members the vendored code does not
// read — polynomial II f*/u* sets, MM/ud/ue — are omitted).
static const mynumber
/**/ d3             = {{0x55555555, IW(0xbfd55555)} }, /* -0.333... */
/**/ d5             = {{IW(0x999997fd), 0x3fc99999} }, /*  0.199... */
/**/ d7             = {{IW(0x923f7603), IW(0xbfc24924)} }, /* -0.142... */
/**/ d9             = {{IW(0xe5129a3b), 0x3fbc71c6} }, /*  0.111... */
/**/ d11            = {{0x22b13c25, IW(0xbfb74580)} }, /* -0.090... */
/**/ d13            = {{IW(0x8b31cbce), 0x3fb375f0} }, /*  0.076... */
/**/ inv16          = {{0x00000000, 0x3fb00000} }, /*  1/16         */
/**/ opi            = {{0x54442d18, 0x400921fb} }, /*  pi           */
/**/ opi1           = {{0x33145c07, 0x3ca1a626} }, /*  pi-opi       */
/**/ mopi           = {{0x54442d18, IW(0xc00921fb)} }, /* -pi           */
/**/ hpi            = {{0x54442d18, 0x3ff921fb} }, /*  pi/2         */
/**/ hpi1           = {{0x33145c07, 0x3c91a626} }, /*  pi/2-hpi     */
/**/ mhpi           = {{0x54442d18, IW(0xbff921fb)} }, /* -pi/2         */
/**/ qpi            = {{0x54442d18, 0x3fe921fb} }, /*  pi/4         */
/**/ mqpi           = {{0x54442d18, IW(0xbfe921fb)} }, /* -pi/4         */
/**/ tqpi           = {{0x7f3321d2, 0x4002d97c} }, /*  3pi/4        */
/**/ mtqpi          = {{0x7f3321d2, IW(0xc002d97c)} }, /* -3pi/4        */
/**/ two500         = {{0x00000000, 0x5f300000} }, /*  2**500       */
/**/ twom500        = {{0x00000000, 0x20b00000} }; /*  2**(-500)    */

#define TWO52 0x1.0p52
#define TWOM1022 0x1.0p-1022

/* Fix the sign and return after stage 1 or stage 2 */
inline double signArctan2(double y, double z) {
  return std::copysign(z, y);
}

double atan2_kernel(double y, double x) {
  int i, de, ux, dx, uy, dy;
  double ax, ay, u, du, v, vv, dv, t1, t2, t3, z, zz, cor;
  mynumber num;

  static const int ep = 59768832,     /*  57*16**5   */
      em = -59768832;                 /* -57*16**5   */

  /* x=NaN or y=NaN */
  num.d = x;
  ux = num.i[HIGH_HALF];
  dx = num.i[LOW_HALF];
  if ((ux & 0x7ff00000) == 0x7ff00000) {
    if (((ux & 0x000fffff) | dx) != 0x00000000)
      return x + y;
  }
  num.d = y;
  uy = num.i[HIGH_HALF];
  dy = num.i[LOW_HALF];
  if ((uy & 0x7ff00000) == 0x7ff00000) {
    if (((uy & 0x000fffff) | dy) != 0x00000000)
      return y + y;
  }

  /* y=+-0 */
  if (uy == 0x00000000) {
    if (dy == 0x00000000) {
      if ((ux & 0x80000000) == 0x00000000)
        return 0;
      else
        return opi.d;
    }
  } else if (uy == IW(0x80000000)) {
    if (dy == 0x00000000) {
      if ((ux & 0x80000000) == 0x00000000)
        return -0.0;
      else
        return mopi.d;
    }
  }

  /* x=+-0 */
  if (x == 0) {
    if ((uy & 0x80000000) == 0x00000000)
      return hpi.d;
    else
      return mhpi.d;
  }

  /* x=+-INF */
  if (ux == 0x7ff00000) {
    if (dx == 0x00000000) {
      if (uy == 0x7ff00000) {
        if (dy == 0x00000000)
          return qpi.d;
      } else if (uy == IW(0xfff00000)) {
        if (dy == 0x00000000)
          return mqpi.d;
      } else {
        if ((uy & 0x80000000) == 0x00000000)
          return 0;
        else
          return -0.0;
      }
    }
  } else if (ux == IW(0xfff00000)) {
    if (dx == 0x00000000) {
      if (uy == 0x7ff00000) {
        if (dy == 0x00000000)
          return tqpi.d;
      } else if (uy == IW(0xfff00000)) {
        if (dy == 0x00000000)
          return mtqpi.d;
      } else {
        if ((uy & 0x80000000) == 0x00000000)
          return opi.d;
        else
          return mopi.d;
      }
    }
  }

  /* y=+-INF */
  if (uy == 0x7ff00000) {
    if (dy == 0x00000000)
      return hpi.d;
  } else if (uy == IW(0xfff00000)) {
    if (dy == 0x00000000)
      return mhpi.d;
  }

  /* either x/y or y/x is very close to zero */
  ax = (x < 0) ? -x : x;
  ay = (y < 0) ? -y : y;
  de = (uy & 0x7ff00000) - (ux & 0x7ff00000);
  if (de >= ep) {
    return ((y > 0) ? hpi.d : mhpi.d);
  } else if (de <= em) {
    if (x > 0) {
      double ret;
      z = ay / ax;
      ret = signArctan2(y, z);
      return ret;
    } else {
      return ((y > 0) ? opi.d : mopi.d);
    }
  }

  /* if either x or y is extremely close to zero, scale abs(x), abs(y). */
  if (ax < twom500.d || ay < twom500.d) {
    ax *= two500.d;
    ay *= two500.d;
  }

  /* Likewise for large x and y.  */
  if (ax > two500.d || ay > two500.d) {
    ax *= twom500.d;
    ay *= twom500.d;
  }

  /* x,y which are neither special nor extreme */
  if (ay < ax) {
    u = ay / ax;
    EMULV(ax, u, v, vv);
    du = ((ay - v) - vv) / ax;
  } else {
    u = ax / ay;
    EMULV(ay, u, v, vv);
    du = ((ax - v) - vv) / ay;
  }

  if (x > 0) {
    /* (i)   x>0, abs(y)< abs(x):  atan(ay/ax) */
    if (ay < ax) {
      if (u < inv16.d) {
        v = u * u;

        zz = du + u * v * (d3.d + v * (d5.d + v * (d7.d + v * (d9.d +
                                                               v * (d11.d + v * d13.d)))));

        z = u + zz;
        /* Max ULP is 0.504.  */
        return signArctan2(y, z);
      }

      i = (TWO52 + 256 * u) - TWO52;
      i -= 16;
      t3 = u - cij[i][0].d;
      EADD(t3, du, v, dv);
      t1 = cij[i][1].d;
      t2 = cij[i][2].d;
      zz = v * t2 +
           (dv * t2 +
            v * v * (cij[i][3].d + v * (cij[i][4].d +
                                        v * (cij[i][5].d + v * cij[i][6].d))));
      z = t1 + zz;
      /* Max ULP is 0.56.  */
      return signArctan2(y, z);
    }

    /* (ii)  x>0, abs(x)<=abs(y):  pi/2-atan(ax/ay) */
    if (u < inv16.d) {
      v = u * u;
      zz = u * v * (d3.d + v * (d5.d + v * (d7.d + v * (d9.d + v * (d11.d +
                                                                     v * d13.d)))));
      ESUB(hpi.d, u, t2, cor);
      t3 = ((hpi1.d + cor) - du) - zz;
      z = t2 + t3;
      /* Max ULP is 0.501.  */
      return signArctan2(y, z);
    }

    i = (TWO52 + 256 * u) - TWO52;
    i -= 16;
    v = (u - cij[i][0].d) + du;

    zz = hpi1.d -
         v * (cij[i][2].d + v * (cij[i][3].d + v * (cij[i][4].d +
                                                    v * (cij[i][5].d + v * cij[i][6].d))));
    t1 = hpi.d - cij[i][1].d;
    z = t1 + zz;
    /* Max ULP is 0.503.  */
    return signArctan2(y, z);
  }

  /* (iii) x<0, abs(x)< abs(y):  pi/2+atan(ax/ay) */
  if (ax < ay) {
    if (u < inv16.d) {
      v = u * u;
      zz = u * v * (d3.d + v * (d5.d + v * (d7.d + v * (d9.d +
                                                        v * (d11.d + v * d13.d)))));
      EADD(hpi.d, u, t2, cor);
      t3 = ((hpi1.d + cor) + du) + zz;
      z = t2 + t3;
      /* Max ULP is 0.501.  */
      return signArctan2(y, z);
    }

    i = (TWO52 + 256 * u) - TWO52;
    i -= 16;
    v = (u - cij[i][0].d) + du;
    zz = hpi1.d +
         v * (cij[i][2].d + v * (cij[i][3].d + v * (cij[i][4].d +
                                                    v * (cij[i][5].d + v * cij[i][6].d))));
    t1 = hpi.d + cij[i][1].d;
    z = t1 + zz;
    /* Max ULP is 0.503.  */
    return signArctan2(y, z);
  }

  /* (iv)  x<0, abs(y)<=abs(x):  pi-atan(ax/ay) */
  if (u < inv16.d) {
    v = u * u;
    zz = u * v * (d3.d + v * (d5.d + v * (d7.d + v * (d9.d + v * (d11.d +
                                                                   v * d13.d)))));
    ESUB(opi.d, u, t2, cor);
    t3 = ((opi1.d + cor) - du) - zz;
    z = t2 + t3;
    /* Max ULP is 0.501.  */
    return signArctan2(y, z);
  }

  i = (TWO52 + 256 * u) - TWO52;
  i -= 16;
  v = (u - cij[i][0].d) + du;
  zz = opi1.d -
       v * (cij[i][2].d + v * (cij[i][3].d + v * (cij[i][4].d +
                                                  v * (cij[i][5].d + v * cij[i][6].d))));
  t1 = opi.d - cij[i][1].d;
  z = t1 + zz;
  /* Max ULP is 0.502.  */
  return signArctan2(y, z);
}

}  // namespace atan2

// ===========================================================================
// asin / acos — e_asin.c (IBM APML uasncs). Only acos is routed by the
// ledger; asin is transcribed for file-level fidelity and proof coverage.
// math_check_force_underflow cannot change returned bits.
// ===========================================================================
namespace asinacos {

// uasncs.h constants, LITTLE_ENDI branch (a1/a2 and c2..c7 are not read by
// the vendored code and are omitted).
static const mynumber
/**/          hp0 = {{0x54442D18, 0x3FF921FB }},  /*  1.5707963267948966     */
/**/          hp1 = {{0x33145C07, 0x3C91A626 }};  /*  6.123233995736766e-17  */

static const double
    f1 = 1.66666666666664110590506577996662E-01,
    f2 = 7.50000000026122686814431784722623E-02,
    f3 = 4.46428561421059750978517350006940E-02,
    f4 = 3.03821268582119319911193410625235E-02,
    f5 = 2.23551211026525610742786300334557E-02,
    f6 = 1.81382903404565056280372531963613E-02;
static const double t24 = 16777216.0, t27 = 134217728.0;
static const double
    rt0 = 9.99999999859990725855365213134618E-01,
    rt1 = 4.99999999495955425917856814202739E-01,
    rt2 = 3.75017500867345182581453026130850E-01,
    rt3 = 3.12523626554518656309172508769531E-01;

double asin_kernel(double x) {
  double x2, xx, res1, p, t, res, r, cor, cc, y, c, z;
  mynumber u, v;
  int4 k, m, n;

  u.x = x;
  m = u.i[HIGH_HALF];
  k = 0x7fffffff & m; /* no sign */

  if (k < 0x3e500000) {
    return x; /* for x->0 => sin(x)=x */
  }
  /*----------------------2^-26 <= |x| < 2^ -3    -----------------*/
  else if (k < 0x3fc00000) {
    x2 = x * x;
    t = ((((((f6 * x2 + f5) * x2 + f4) * x2 + f3) * x2 + f2) * x2 + f1) *
         (x2 * x));
    res = x + t; /*  res=arcsin(x) according to Taylor series  */
    return res;
  }
  /*---------------------0.125 <= |x| < 0.5 -----------------------------*/
  else if (k < 0x3fe00000) {
    if (k < 0x3fd00000)
      n = 11 * ((k & 0x000fffff) >> 15);
    else
      n = 11 * ((k & 0x000fffff) >> 14) + 352;
    if (m > 0)
      xx = x - asncs.x[n];
    else
      xx = -x - asncs.x[n];
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] + xx * asncs.x[n + 6])))) +
        asncs.x[n + 7];
    t += p;
    res = asncs.x[n + 8] + t;
    return (m > 0) ? res : -res;
  } /*   else  if (k < 0x3fe00000)    */
  /*-------------------- 0.5 <= |x| < 0.75 -----------------------------*/
  else if (k < 0x3fe80000) {
    n = 1056 + ((k & 0x000fe000) >> 11) * 3;
    if (m > 0)
      xx = x - asncs.x[n];
    else
      xx = -x - asncs.x[n];
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * asncs.x[n + 7]))))) +
        asncs.x[n + 8];
    t += p;
    res = asncs.x[n + 9] + t;
    return (m > 0) ? res : -res;
  } /*   else  if (k < 0x3fe80000)    */
  /*--------------------- 0.75 <= |x|< 0.921875 ----------------------*/
  else if (k < 0x3fed8000) {
    n = 992 + ((k & 0x000fe000) >> 13) * 13;
    if (m > 0)
      xx = x - asncs.x[n];
    else
      xx = -x - asncs.x[n];
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * (asncs.x[n + 7] +
                                                 xx * asncs.x[n + 8])))))) +
        asncs.x[n + 9];
    t += p;
    res = asncs.x[n + 10] + t;
    return (m > 0) ? res : -res;
  } /*   else  if (k < 0x3fed8000)    */
  /*-------------------0.921875 <= |x| < 0.953125 ------------------------*/
  else if (k < 0x3fee8000) {
    n = 884 + ((k & 0x000fe000) >> 13) * 14;
    if (m > 0)
      xx = x - asncs.x[n];
    else
      xx = -x - asncs.x[n];
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * (asncs.x[n + 7] +
                                                 xx * (asncs.x[n + 8] +
                                                       xx * asncs.x[n + 9]))))))) +
        asncs.x[n + 10];
    t += p;
    res = asncs.x[n + 11] + t;
    return (m > 0) ? res : -res;
  }

  /*--------------------0.953125 <= |x| < 0.96875 ------------------------*/
  else if (k < 0x3fef0000) {
    n = 768 + ((k & 0x000fe000) >> 13) * 15;
    if (m > 0)
      xx = x - asncs.x[n];
    else
      xx = -x - asncs.x[n];
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * (asncs.x[n + 7] +
                                                 xx * (asncs.x[n + 8] +
                                                       xx * (asncs.x[n + 9] +
                                                             xx * asncs.x[n + 10])))))))) +
        asncs.x[n + 11];
    t += p;
    res = asncs.x[n + 12] + t;
    return (m > 0) ? res : -res;
  }
  /*--------------------0.96875 <= |x| < 1 --------------------------------*/
  else if (k < 0x3ff00000) {
    z = 0.5 * ((m > 0) ? (1.0 - x) : (1.0 + x));
    v.x = z;
    k = v.i[HIGH_HALF];
    t = inroot[(k & 0x001fffff) >> 14] * powtwo[511 - (k >> 21)];
    r = 1.0 - t * t * z;
    t = t * (rt0 + r * (rt1 + r * (rt2 + r * rt3)));
    c = t * z;
    t = c * (1.5 - 0.5 * t * c);
    y = (c + t24) - t24;
    cc = (z - y * y) / (t + y);
    p = ((((((f6 * z + f5) * z + f4) * z + f3) * z + f2) * z + f1) * z);
    cor = (hp1.x - 2.0 * cc) - 2.0 * (y + cc) * p;
    res1 = hp0.x - 2.0 * y;
    res = res1 + cor;
    return (m > 0) ? res : -res;
  } /*   else  if (k < 0x3ff00000)    */
  /*---------------------------- |x|>=1 -------------------------------*/
  else if (k == 0x3ff00000 && u.i[LOW_HALF] == 0)
    return (m > 0) ? hp0.x : -hp0.x;
  else
    return (x - x) / (x - x);
}

double acos_kernel(double x) {
  double x2, xx, res1, p, t, res, r, cor, cc, y, c, z;
  mynumber u, v;
  int4 k, m, n;
  u.x = x;
  m = u.i[HIGH_HALF];
  k = 0x7fffffff & m;
  /*-------------------  |x|<2.77556*10^-17 ----------------------*/
  if (k < 0x3c880000)
    return hp0.x;

  /*-----------------  2.77556*10^-17 <= |x| < 2^-3 --------------*/
  else if (k < 0x3fc00000) {
    x2 = x * x;
    t = ((((((f6 * x2 + f5) * x2 + f4) * x2 + f3) * x2 + f2) * x2 + f1) *
         (x2 * x));
    r = hp0.x - x;
    cor = (((hp0.x - r) - x) + hp1.x) - t;
    res = r + cor;
    return res;
  } /*   else  if (k < 0x3fc00000)    */
  /*----------------------  0.125 <= |x| < 0.5 --------------------*/
  else if (k < 0x3fe00000) {
    if (k < 0x3fd00000)
      n = 11 * ((k & 0x000fffff) >> 15);
    else
      n = 11 * ((k & 0x000fffff) >> 14) + 352;
    if (m > 0)
      xx = x - asncs.x[n];
    else
      xx = -x - asncs.x[n];
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] + xx * asncs.x[n + 6])))) +
        asncs.x[n + 7];
    t += p;
    y = (m > 0) ? (hp0.x - asncs.x[n + 8]) : (hp0.x + asncs.x[n + 8]);
    t = (m > 0) ? (hp1.x - t) : (hp1.x + t);
    res = y + t;
    return res;
  }

  /*--------------------------- 0.5 <= |x| < 0.75 ---------------------*/
  else if (k < 0x3fe80000) {
    n = 1056 + ((k & 0x000fe000) >> 11) * 3;
    if (m > 0) {
      xx = x - asncs.x[n];
    } else {
      xx = -x - asncs.x[n];
    }
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * asncs.x[n + 7]))))) +
        asncs.x[n + 8];
    t += p;
    y = (m > 0) ? (hp0.x - asncs.x[n + 9]) : (hp0.x + asncs.x[n + 9]);
    t = (m > 0) ? (hp1.x - t) : (hp1.x + t);
    res = y + t;
    return res;
  }

  /*------------------------- 0.75 <= |x| < 0.921875 -------------*/
  else if (k < 0x3fed8000) {
    n = 992 + ((k & 0x000fe000) >> 13) * 13;
    if (m > 0) {
      xx = x - asncs.x[n];
    } else {
      xx = -x - asncs.x[n];
    }
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * (asncs.x[n + 7] +
                                                 xx * asncs.x[n + 8])))))) +
        asncs.x[n + 9];
    t += p;
    y = (m > 0) ? (hp0.x - asncs.x[n + 10]) : (hp0.x + asncs.x[n + 10]);
    t = (m > 0) ? (hp1.x - t) : (hp1.x + t);
    res = y + t;
    return res;
  }

  /*-------------------0.921875 <= |x| < 0.953125 ------------------*/
  else if (k < 0x3fee8000) {
    n = 884 + ((k & 0x000fe000) >> 13) * 14;
    if (m > 0) {
      xx = x - asncs.x[n];
    } else {
      xx = -x - asncs.x[n];
    }
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * (asncs.x[n + 7] +
                                                 xx * (asncs.x[n + 8] +
                                                       xx * asncs.x[n + 9]))))))) +
        asncs.x[n + 10];
    t += p;
    y = (m > 0) ? (hp0.x - asncs.x[n + 11]) : (hp0.x + asncs.x[n + 11]);
    t = (m > 0) ? (hp1.x - t) : (hp1.x + t);
    res = y + t;
    return res;
  }

  /*--------------------0.953125 <= |x| < 0.96875 ----------------*/
  else if (k < 0x3fef0000) {
    n = 768 + ((k & 0x000fe000) >> 13) * 15;
    if (m > 0) {
      xx = x - asncs.x[n];
    } else {
      xx = -x - asncs.x[n];
    }
    t = asncs.x[n + 1] * xx;
    p = xx * xx * (asncs.x[n + 2] +
                   xx * (asncs.x[n + 3] +
                         xx * (asncs.x[n + 4] +
                               xx * (asncs.x[n + 5] +
                                     xx * (asncs.x[n + 6] +
                                           xx * (asncs.x[n + 7] +
                                                 xx * (asncs.x[n + 8] +
                                                       xx * (asncs.x[n + 9] +
                                                             xx * asncs.x[n + 10])))))))) +
        asncs.x[n + 11];
    t += p;
    y = (m > 0) ? (hp0.x - asncs.x[n + 12]) : (hp0.x + asncs.x[n + 12]);
    t = (m > 0) ? (hp1.x - t) : (hp1.x + t);
    res = y + t;
    return res;
  }
  /*-----------------0.96875 <= |x| < 1 ---------------------------*/

  else if (k < 0x3ff00000) {
    z = 0.5 * ((m > 0) ? (1.0 - x) : (1.0 + x));
    v.x = z;
    k = v.i[HIGH_HALF];
    t = inroot[(k & 0x001fffff) >> 14] * powtwo[511 - (k >> 21)];
    r = 1.0 - t * t * z;
    t = t * (rt0 + r * (rt1 + r * (rt2 + r * rt3)));
    c = t * z;
    t = c * (1.5 - 0.5 * t * c);
    y = (t27 * c + c) - t27 * c;
    cc = (z - y * y) / (t + y);
    p = ((((((f6 * z + f5) * z + f4) * z + f3) * z + f2) * z + f1) * z);
    if (m < 0) {
      cor = (hp1.x - cc) - (y + cc) * p;
      res1 = hp0.x - y;
      res = res1 + cor;
      return (res + res);
    } else {
      cor = cc + p * (y + cc);
      res = y + cor;
      return (res + res);
    }
  } /*   else  if (k < 0x3ff00000)    */

  /*---------------------------- |x|>=1 -----------------------*/
  else if (k == 0x3ff00000 && u.i[LOW_HALF] == 0)
    return (m > 0) ? 0 : 2.0 * hp0.x;
  else
    return (x - x) / (x - x);
}

}  // namespace asinacos

}  // anonymous namespace

}  // namespace kalos::detail::hostlineage_src

// ---------------------------------------------------------------------------
// Public API (kalos::detail) — thin dispatchers into the vendored kernels.
// ---------------------------------------------------------------------------
namespace kalos::detail {

double glibc_sin(double x) { return hostlineage_src::sincos::sin_kernel(x); }
double glibc_cos(double x) { return hostlineage_src::sincos::cos_kernel(x); }
double glibc_exp(double x) { return hostlineage_src::exp::exp_kernel(x); }
double glibc_log(double x) { return hostlineage_src::log::log_kernel(x); }
double glibc_pow(double x, double y) {
  return hostlineage_src::pow::pow_kernel(x, y);
}
double glibc_hypot(double x, double y) {
  return hostlineage_src::hypotk::hypot_kernel(x, y);
}
double glibc_atan2(double y, double x) {
  return hostlineage_src::atan2::atan2_kernel(y, x);
}
double glibc_acos(double x) {
  return hostlineage_src::asinacos::acos_kernel(x);
}
double glibc_asin(double x) {
  return hostlineage_src::asinacos::asin_kernel(x);
}

}  // namespace kalos::detail
