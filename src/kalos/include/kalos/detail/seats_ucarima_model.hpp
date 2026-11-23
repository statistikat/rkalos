#pragma once

// kalos M5 slice 24-C seam: the structured UCARIMA model mirror shared by
// 24-B (decompose_ucarima.cpp fact emission) and 24-C (the WK extraction
// machine). This is the EXACT code 24-B carried in its anonymous namespace
// (byte-for-byte, with `inline' added for header linkage), moved verbatim so
// the two slices consume one and the same pinned mirror of
// jdplus.toolkit.base.core.arima.ArimaModel +
// jdplus.toolkit.base.core.ucarima.UcarimaModel: the SymmetricFilter
// algebra, UnitRootsSolver, the EV2/EV1 factorizer chain
// (SymmetricFilter.factorize), the BackFilter/Polynomial simplifying tools,
// the ArimaModel spectral factorization cache (getMa/getIV) and the
// setVarianceMax/compactUcarimaModel pipeline. Behavior is unchanged from
// 24-B by construction; the 24-B corpus goldens are the regression gate.
//
// Names resolve exactly as they did inside `namespace kalos { namespace {'
// of decompose_ucarima.cpp: from kalos::detail::seats_model the outward
// lookup finds kalos::detail's engine helpers (Poly, SymF, cplx, the *_g /
// *_l comparison shims) exactly as before.

#include <exception>
#include "kalos/detail/hostlineage.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "kalos/detail/ucarima_engine.hpp"

namespace kalos {
namespace detail {
namespace seats_model {

using detail::Poly;
using detail::SymF;
using detail::cplx;
using vd = std::vector<double>;
using vi = std::vector<int>;

using detail::c_abs;
using detail::c_cart;
using detail::c_conj;
using detail::c_inv;
using detail::c_minus;
using detail::c_scale;
using detail::c_times;
using detail::ge_g;
using detail::ge_l;
using detail::gt_g;
using detail::gt_l;
using detail::is_nan;
using detail::kFeEps;
using detail::kMax;
using detail::kNan;
using detail::kPi;
using detail::leja_order;
using detail::lt_g;
using detail::lt_l;
using detail::p_from_roots;
using detail::p_minus_p;
using detail::p_of;
using detail::p_quotient;
using detail::p_raw;
using detail::p_times;
using detail::p_times_d;
using detail::pl_divide;
using detail::pl_times;
using detail::sym_convolution;
using detail::trim_copy;

inline std::size_t idx(int i) { return static_cast<std::size_t>(i); }
inline std::size_t idx(std::size_t i) { return i; }

// The probe's g(v) == String.format(Locale.ROOT, "%.17g", v). Java's %g
// differs from C's in two ways that matter here: it NEVER strips trailing
// zeros (it pads the SHORTEST round-trip digits - Double.toString's - up to
// `precision' significant digits), and its scientific switch is
// exp < -4 || exp >= precision. Reproduced on top of std::to_chars'
// shortest round-trip digits (the same canonical digit string Java's
// Dragonbox-based Double.toString produces).
inline std::string fmt_g(double v) {
  if (std::isnan(v)) return "NaN";
  if (std::isinf(v)) return v < 0 ? "-Infinity" : "Infinity";
  char buf[64];
  auto res = std::to_chars(buf, buf + sizeof buf, v, std::chars_format::scientific);
  std::string s(buf, static_cast<std::size_t>(res.ptr - buf));
  std::string out;
  if (!s.empty() && s[0] == '-') {
    out = "-";
    s.erase(0, 1);
  }
  std::string mant = s.substr(0, s.find('e'));
  int exp10 = 0;
  if (mant.size() != s.size()) exp10 = std::stoi(s.substr(mant.size() + 1));
  std::string digits;
  for (char c : mant)
    if (c != '.') digits += c;
  if (digits.empty() || digits.find_first_not_of('0') == std::string::npos) {
    digits = "0";
    exp10 = 0;
  }
  constexpr int p = 17;
  if (exp10 < -4 || exp10 >= p) {  // scientific
    out += digits[0];
    out += '.';
    std::string frac = digits.substr(1);
    frac.append(static_cast<std::size_t>(p - 1) - frac.size(), '0');
    out += frac;
    out += 'e';
    out += exp10 < 0 ? '-' : '+';
    char eb[16];
    std::snprintf(eb, sizeof eb, "%d", exp10 < 0 ? -exp10 : exp10);
    if (std::strlen(eb) < 2) out += '0';
    out += eb;
  } else {  // fixed, p significant digits (shortest is <= p, so pad only)
    digits.append(static_cast<std::size_t>(p) - digits.size(), '0');
    if (exp10 < 0) {
      out += "0.";
      out.append(static_cast<std::size_t>(-exp10 - 1), '0');
      out += digits;
    } else {
      std::size_t ip = static_cast<std::size_t>(exp10) + 1;
      out += digits.substr(0, ip);
      if (digits.size() > ip) {
        out += '.';
        out += digits.substr(ip);
      }
    }
  }
  return out;
}

inline std::string coeffs_row(const vd& c) {
  std::string s;
  for (std::size_t i = 0; i < c.size(); ++i) {
    if (i) s += ' ';
    s += fmt_g(c[i]);
  }
  return s;
}

// ---------------------------------------------------------------------------
// Complex helpers the engine does not need (ComplexType component ops).
// ---------------------------------------------------------------------------
inline cplx c_add(cplx a, cplx b) { return {a.re + b.re, a.im + b.im}; }
inline cplx c_real(double x) { return {x, 0.0}; }

// Support.sign(p, r) (jdplus matrices.Support, pinned bytecode):
//   r >= 0 or NaN: p < 0 ? -p : p    (== |p|, NaN -> p)
//   r < 0        : p < 0 ?  p : -p   (== -|p|)
// NOT copysign: the comparison flavor is plain dcmpl and -0.0 counts as r >= 0.
inline double sign_like(double p, double r) {
  if (r < 0.0) return (p < 0.0) ? p : -p;
  return (p < 0.0) ? -p : p;
}

// ComplexMath.sqrt(x + iy) (ComplexMath.txt @1045-1090). Branch tests are
// dcmpl: NaN r falls THROUGH to the re test (no ZERO), NaN re takes the
// re > 0 branch; the im-flip is dcmpg (NaN skips it). The re <= 0 branch
// ends in cart((0.5*im)/t, t) -- REAL part first is (0.5*im)/t (bytecode
// @71-81; an earlier decode had cart(t, ...) and produced non-square-roots
// -- corrected 2026-11-21).
inline cplx c_sqrt(cplx z) {
  double r = c_abs(z);
  if (r <= 0.0) return {0.0, 0.0};
  if (!(z.re <= 0.0)) {  // dcmpl ifle: NaN re lands here
    double t = std::sqrt(0.5 * (r + z.re));
    return {t, (0.5 * z.im) / t};
  }
  double t = std::sqrt(0.5 * (r - z.re));
  if (z.im < 0.0) t = -t;  // dcmpg ifge
  return {(0.5 * z.im) / t, t};
}

// JVM exceptions modelled for control flow (caught exactly where the mirrored
// class catches them; anywhere else they surface at the probe's try/catch
// sites, which print the pinned class name).
struct MaError : std::exception {};    // ArimaException("invalid ARIMA")
struct JvmNpe : std::exception {};     // java.lang.NullPointerException
struct Aioobe : std::exception {};     // java.lang.ArrayIndexOutOfBounds

inline const char* exception_class(const std::exception_ptr& e) {
  try {
    std::rethrow_exception(e);
  } catch (const MaError&) {
    return "jdplus.toolkit.base.core.arima.ArimaException";
  } catch (const JvmNpe&) {
    return "java.lang.NullPointerException";
  } catch (const Aioobe&) {
    return "java.lang.ArrayIndexOutOfBoundsException";
  } catch (const detail::jvm_decompose_abort&) {
    return "jdplus.toolkit.base.core.math.matrices.MatrixException";
  } catch (const std::overflow_error&) {
    return "java.lang.ArithmeticException";
  } catch (...) {
    return "java.lang.RuntimeException";
  }
}

// ---------------------------------------------------------------------------
// Polynomial.plus(Polynomial) (q zero -> this; else Polynomials.plus +
// Coefficients.ofInternal trim).
// ---------------------------------------------------------------------------
vd pl_plus(const vd& a, const vd& b) {  // Polynomials.plus
  vd r(std::max(a.size(), b.size()), 0.0);
  for (std::size_t i = 0; i < a.size(); ++i) r[i] += a[i];
  for (std::size_t i = 0; i < b.size(); ++i) r[i] += b[i];
  return r;
}

Poly p_plus(const Poly& a, const Poly& b) {  // Polynomial.plus
  if (b.is_zero()) return a;
  return p_raw(trim_copy(pl_plus(a.c, b.c)));
}

// Polynomial scalar ops on raw coefficient arrays (SymmetricFilter delegates
// to exactly these; the results reach `new SymmetricFilter(Polynomial)`
// UNLESS the delegate trims, which the comments spell out).
vd pl_times_d(const vd& a, double d) {  // Polynomial.times(double)
  if (d == 0.0 || (a.size() == 1 && a[0] == 0.0)) return {0.0};
  if (d == 1.0) return a;
  vd r = a;
  for (double& x : r) x *= d;
  return r;
}

bool sf_is_null(const SymF& f) { return f.a.size() == 1 && f.a[0] == 0.0; }

// SymmetricFilter instance arithmetic (every method is "new
// SymmetricFilter(polynomial.op(...))" -- the raw constructor, no trim,
// except times(SF), which ends in ofInternal (trim) as pinned below).
SymF sf_times_d(const SymF& f, double d) { return SymF(pl_times_d(f.a, d)); }
SymF sf_multiply_d(double d, const SymF& f) { return sf_times_d(f, d); }

SymF sf_plus(const SymF& a, const SymF& b) {  // plus(SF): trims
  return SymF(p_plus(Poly(a.a), Poly(b.a)).c);
}

SymF sf_minus(const SymF& a, const SymF& b) {  // minus(SF): trims iff equal
  return SymF(p_minus_p(Poly(a.a), Poly(b.a)).c);  // degrees
}

// SymmetricFilter.times(SymmetricFilter) (pinned byte-exact, including the
// cross-term folding with the 2.0 double-count and the terminating
// ofInternal (trim)).
inline SymF sf_times(const SymF& x, const SymF& y) {
  int n = static_cast<int>(x.a.size()) - 1;
  int m = static_cast<int>(y.a.size()) - 1;
  const vd& pa = x.a;
  const vd& pb = y.a;
  vd r(static_cast<std::size_t>(n + m + 1), 0.0);
  if (pb[0] != 0.0)  // dcmpl ifeq: NaN keeps the loop
    for (int i = 0; i <= n; ++i) r[idx(i)] += pa[idx(i)] * pb[0];
  if (pa[0] != 0.0)
    for (int i = 1; i <= m; ++i) r[idx(i)] += pb[idx(i)] * pa[0];
  for (int i = 1; i <= n; ++i) {
    if (pa[idx(i)] == 0.0) continue;
    for (int j = 1; j <= m; ++j) {
      if (pb[idx(j)] == 0.0) continue;
      double w = pa[idx(i)] * pb[idx(j)];
      r[idx(i + j)] += w;
      if (i > j)
        r[idx(i - j)] += w;
      else if (i < j)
        r[idx(j - i)] += w;
      else
        r[0] += 2.0 * w;
    }
  }
  return SymF(trim_copy(r));  // ofInternal
}

// SymmetricFilter.weights(): t -> a[|t|]; weightsToArray(): the palindrome
// [-deg .. deg].
inline double sf_weight(const SymF& f, int t) {
  return f.a[idx(t < 0 ? -t : t)];
}

inline vd sf_weights_array(const SymF& f) {
  int d = static_cast<int>(f.a.size()) - 1;
  vd w(static_cast<std::size_t>(2 * d + 1));
  for (int t = -d; t <= d; ++t) w[idx(t + d)] = sf_weight(f, t);
  return w;
}

// ---------------------------------------------------------------------------
// jdplus ...core.utils.IVector (Java capacity semantics): add() grows to
// sz+4, add(int[],m) reallocs to ((sz+m)/4+1)*4 slots, clear() drops the
// array, compact() KEEPS the capacity (the entries beyond sz keep their
// previous values - IVector.sqrt reads one of them on odd sizes, which the
// mirror reproduces literally).
// ---------------------------------------------------------------------------
struct IVec {
  vi vals;
  int sz = 0;

  int size() const { return sz; }
  bool is_empty() const { return sz == 0; }

  void add(int v) {  // IVector.add(int)
    if (vals.empty()) {
      vals.assign(4, 0);
      vals[0] = v;
      sz = 1;
      return;
    }
    if (static_cast<int>(vals.size()) == sz) {
      vi n(static_cast<std::size_t>(sz) + 4, 0);
      std::copy(vals.begin(), vals.begin() + sz, n.begin());
      n[idx(sz)] = v;
      ++sz;
      vals = std::move(n);
      return;
    }
    vals[idx(sz++)] = v;
  }
  void add(const vi& a, int m) {  // IVector.add(int[], int)
    if (!vals.empty() && sz + m <= static_cast<int>(vals.size())) {
      for (int i = 0; i < m; ++i) vals[idx(sz++)] = a[idx(i)];
      return;
    }
    int nsz = ((sz + m) / 4 + 1) * 4;
    vi n(static_cast<std::size_t>(nsz), 0);
    for (int i = 0; i < sz; ++i) n[idx(i)] = vals[idx(i)];
    for (int i = 0; i < m; ++i) n[idx(sz++)] = a[idx(i)];
    vals = std::move(n);
  }
  void add(const IVec& o) {
    if (o.sz != 0) add(o.vals, o.sz);
  }
  void clear() {
    sz = 0;
    vals.clear();
  }
  // IVector.compact(): drop zero entries inside [0, sz); the capacity (and
  // the stale values beyond the retained entries) is kept.
  void compact() {
    int k = 0;
    for (int i = 0; i < sz; ++i)
      if (vals[idx(i)] != 0) vals[idx(k++)] = vals[idx(i)];
    sz = k;
  }
  void rsort_desc() {  // IVector.RSort (descending)
    std::sort(vals.begin(), vals.begin() + sz, std::greater<int>());
  }

  static IVec from(const vi& a) {
    IVec v;
    v.vals = a;
    v.sz = static_cast<int>(a.size());
    return v;
  }
};

// IVector.simplify(v0, v1) (pinned): RSort both descending, two-pointer
// cross-cancellation (equal -> zero both and advance both; else v0[j] >
// v1[i] -> ++j, else ++i); compact both on any hit.
inline void iv_simplify(IVec& v0, IVec& v1) {
  v0.rsort_desc();
  v1.rsort_desc();
  int i = 0, j = 0, cnt = 0;
  while (i < v1.sz && j < v0.sz) {
    if (v1.vals[idx(i)] == v0.vals[idx(j)]) {
      ++cnt;
      v1.vals[idx(i++)] = 0;
      v0.vals[idx(j++)] = 0;
    } else if (v0.vals[idx(j)] > v1.vals[idx(i)]) {
      ++j;
    } else {
      ++i;
    }
  }
  if (cnt > 0) {
    v0.compact();
    v1.compact();
  }
}

// IVector.sqrt(): pairwise equality; null on mismatch; an odd size reads the
// stale capacity slot (or throws AIOOBE when size == capacity). The JVM
// AIOOBE is signalled here so EigenValuesDecomposer2.decompose() catches it.
struct SqrtResult {
  bool threw = false;
  bool ok = false;
  IVec out;
};

inline SqrtResult iv_sqrt(const IVec& v) {
  SqrtResult r;
  IVec res;
  for (int i = 0; i < v.sz; i += 2) {
    if (i + 1 >= static_cast<int>(v.vals.size())) {
      r.threw = true;  // ArrayIndexOutOfBoundsException
      return r;
    }
    if (v.vals[idx(i)] != v.vals[idx(i + 1)]) return r;  // null (mismatch)
    res.add(v.vals[idx(i)]);
  }
  r.ok = true;
  r.out = std::move(res);
  return r;
}

// ---------------------------------------------------------------------------
// jdplus ...core.math.polynomials.UnitRoots + UnitRootsSolver.
// m_n/m_d are IVector LISTS of the (1 - B^m) factor exponents.
// ---------------------------------------------------------------------------
void ur_int_mul(vi& a, int& len, int m) {  // UnitRoots.mul(int[], int, int)
  for (int t = len - 1; t >= 0; --t) a[idx(t + m)] -= a[idx(t)];
  len += m;
}

int ur_int_div(vi& a, int len, int m) {  // UnitRoots.div(int[], int, int)
  for (int t = len - 1; t >= m; --t) {
    a[idx(t - m)] += a[idx(t)];
    a[idx(t)] = -a[idx(t)];  // verbatim JVM quirk
  }
  for (int t = 0; t < m && t < len; ++t)
    if (a[idx(t)] != 0) return -1;
  for (int t = m; t < len; ++t) a[idx(t - m)] = a[idx(t)];
  return len - m;
}

void ur_dbl_mul(vd& a, int& len, int m) {  // UnitRoots.mul(double[], int, int)
  for (int t = len - 1; t >= 0; --t) a[idx(t + m)] -= a[idx(t)];
  len += m;
}

int ur_dbl_div(vd& a, int len, int m) {  // UnitRoots.div(double[], int, int)
  for (int t = m; t < len; ++t) a[idx(t)] += a[idx(t - m)];
  for (int t = len - m; t < len; ++t)
    if (std::fabs(a[idx(t)]) > 1.0E-8) return -1;  // NaN tolerated (dcmpg)
  return len - m;
}

struct UnitRoots {
  IVec m_n, m_d;

  int roots_count() const {  // getRootsCount()
    int s = 0;
    for (int i = 0; i < m_n.sz; ++i) s += m_n.vals[idx(i)];
    for (int i = 0; i < m_d.sz; ++i) s -= m_d.vals[idx(i)];
    return s;
  }
  int denom_degree() const {
    int s = 0;
    for (int i = 0; i < m_d.sz; ++i) s += m_d.vals[idx(i)];
    return s;
  }
  int num_degree() const {
    int s = 0;
    for (int i = 0; i < m_n.sz; ++i) s += m_n.vals[idx(i)];
    return s;
  }

  void add(int m) { m_n.add(m); }  // UnitRoots.add(int)
  void add(const UnitRoots& o) {   // UnitRoots.add(UnitRoots)
    m_n.add(o.m_n);
    m_d.add(o.m_d);
  }
  void clear() {
    m_n.clear();
    m_d.clear();
  }

  // UnitRoots.addOnly(int[], int): push divs[i] into m_n, and for every
  // smaller divs[j] dividing it, recursively push divs[j] into m_d.
  void addonly(const vi& divs, int i) {
    int v = divs[idx(i)];
    m_n.add(v);
    for (int j = i - 1; j >= 0; --j)
      if (v % divs[idx(j)] == 0) ddelonly(divs, j);
  }
  void ddelonly(const vi& divs, int j) {
    int v = divs[idx(j)];
    m_d.add(v);
    for (int k = j - 1; k >= 0; --k)
      if (v % divs[idx(k)] == 0) addonly(divs, k);
  }
  // UnitRoots.removeOnly(int[], int)
  void remove_only(const vi& divs, int i) {
    int v = divs[idx(i)];
    m_d.add(v);
    for (int j = i - 1; j >= 0; --j)
      if (v % divs[idx(j)] == 0) addonly(divs, j);
  }

  // UnitRoots.simplify() + IVector.sqrt() pair-collapse (the receiver copy
  // is simplified; nothing downstream re-reads the solver's roots).
  std::optional<UnitRoots> sqrt_() const {
    UnitRoots u = *this;
    iv_simplify(u.m_n, u.m_d);
    SqrtResult n = iv_sqrt(u.m_n);
    if (n.threw) throw Aioobe{};
    if (!n.ok) return std::nullopt;
    SqrtResult d = iv_sqrt(u.m_d);
    if (d.threw) throw Aioobe{};
    if (!d.ok) return std::nullopt;
    UnitRoots r;
    r.m_n = std::move(n.out);
    r.m_d = std::move(d.out);
    return r;
  }

  // UnitRoots.coefficients() + asPolynomial(): exact integer filter algebra.
  // A negative length from div() is the JVM's AIOOBE (mirrored).
  Poly as_polynomial() const {
    vi a(static_cast<std::size_t>(num_degree()) + 1, 0);
    int len = 1;
    a[0] = 1;
    for (int i = 0; i < m_n.sz; ++i) ur_int_mul(a, len, m_n.vals[idx(i)]);
    for (int i = 0; i < m_d.sz; ++i) {
      len = ur_int_div(a, len, m_d.vals[idx(i)]);
      if (len < 0) throw Aioobe{};
    }
    vd out(static_cast<std::size_t>(len));
    for (int i = 0; i < len; ++i) out[idx(i)] = static_cast<double>(a[idx(i)]);
    return p_of(std::move(out));
  }

  // UnitRoots.divide(Polynomial, UnitRoots): null on inexact division.
  static const Poly* divide(const Poly& p, const UnitRoots& u, Poly& store) {
    if (p.degree() < u.roots_count()) return nullptr;
    vd a(static_cast<std::size_t>(p.degree() + 1 + u.denom_degree()), 0.0);
    for (std::size_t i = 0; i < p.c.size(); ++i) a[i] = p.c[i];
    int len = p.degree() + 1;
    for (int i = 0; i < u.m_d.sz; ++i) ur_dbl_mul(a, len, u.m_d.vals[idx(i)]);
    for (int i = 0; i < u.m_n.sz; ++i) {
      len = ur_dbl_div(a, len, u.m_n.vals[idx(i)]);
      if (len < 0) return nullptr;
    }
    a.resize(static_cast<std::size_t>(len));
    store = p_of(std::move(a));
    return &store;
  }
};

// Arithmetics.divisors(n, out): 0 when n == 1; else every d >= 1 with
// 2d <= n and n % d == 0, ascending.
inline int ur_divisors(int n, vi& out) {
  out.clear();
  if (n == 1) return 0;
  for (int d = 1; 2 * d <= n; ++d)
    if (n % d == 0) out.push_back(d);
  return static_cast<int>(out.size());
}

// UnitRootsSolver(0) state after factorize: the roots plus the solver's
// remainder() (the quotient chain result q). 24-C.1 cluster 3 addition:
// ArimaModel.simplifyUr reads solver.remainder(), so the loop is exposed
// whole; ur_solver_factorize keeps its original roots-only contract.
struct UrSolverResult {
  UnitRoots roots;
  Poly remainder;
};

inline UrSolverResult ur_solver_run(const Poly& p) {
  UnitRoots ur;
  Poly q = p;
  int n = 2 * p.degree() + 1;
  vi divs;
  UnitRoots tmp;
  auto pnorm = [&q](int nn) {
    // Raw libm for the cos/sin below (CI-FDLIBM ledger 2026-11-21): the dcmpg
    // port's Math.cos/sin reaches the JVM's live C2 dcos/dsin stubs
    // (IBM-APML lineage == host libm), and the differencing decisions of every
    // qualified platform (incl. macOS arm64) replay bit-exactly on host-libm
    // bits; the FdLibm bytecode kernels would move the pnorm boundary. See
    // include/kalos/detail/ucarima_engine.hpp (c_arg note) for the probe data.
    if (nn == 1) return std::fabs(detail::p_value_at(q, 1.0));
    if (nn == 2) return std::fabs(detail::p_value_at(q, -1.0));
    double th = 6.283185307179586 / static_cast<double>(nn);
    cplx z{kalos::detail::glibc_cos(th), kalos::detail::glibc_sin(th)};
    cplx acc{0.0, 0.0};
    for (int i = q.degree(); i >= 0; --i) acc = c_add(c_times(acc, z), c_real(q.c[idx(i)]));
    return c_abs(acc);
  };
  while (n >= 1) {
    if (!(pnorm(n) < 1.0E-9)) {  // dcmpg
      --n;
      continue;
    }
    tmp.add(n);
    int cnt = ur_divisors(n, divs);
    for (int i = cnt - 1; i >= 0; --i)
      if (ge_g(pnorm(divs[idx(i)]), 1.0E-9))  // dcmpl ifle: NaN removes
        tmp.remove_only(divs, i);
    Poly store;
    const Poly* w = UnitRoots::divide(q, tmp, store);
    if (w != nullptr) {
      q = *w;
      ur.add(tmp);
    } else {
      --n;
    }
    tmp.clear();
  }
  return {std::move(ur), std::move(q)};
}

inline UnitRoots ur_solver_factorize(const Poly& p) {
  return ur_solver_run(p).roots;
}

// ---------------------------------------------------------------------------
// The eigen pipeline (EV2's roots(double[])): companion matrix ->
// GeneralEigenSystem (balance; companion matrices are Hessenberg, the
// isHessenberg check is ported anyway) -> EigenRoutines.hessenbergQR, all
// column-major (element (r, c) at r + c*n), verbatim from the decoded
// dumps (tmp_opencode/hqr-decoded.md).
// ---------------------------------------------------------------------------
struct CMat {
  int n;
  vd a;  // a[i + j*n]
  explicit CMat(int n_) : n(n_), a(static_cast<std::size_t>(n_) * static_cast<std::size_t>(n_), 0.0) {}
  double at(int i, int j) const { return a[idx(i) + idx(j) * idx(n)]; }
  double& ref(int i, int j) { return a[idx(i) + idx(j) * idx(n)]; }
};

// EigenRoutines.balance (EISPACK balanc, m_radix = 2). Pinned: the contiguous
// block a[i*n + k] (COLUMN i) is csum; the strided a[i + k*n] (ROW i) is
// rsum. Row scaling (1/s) is applied BEFORE the column scaling (s).
inline void bal(CMat& m) {
  int n = m.n;
  constexpr double radix2 = 4.0;
  bool again = true;
  while (again) {
    again = false;
    for (int i = 0; i < n; ++i) {
      double csum = 0.0, rsum = 0.0;
      for (int k = 0; k < n; ++k) {
        if (k == i) continue;
        csum += std::fabs(m.a[idx(i) * idx(n) + idx(k)]);  // column i
        rsum += std::fabs(m.a[idx(i) + idx(k) * idx(n)]);  // row i
      }
      if (rsum == 0.0 || csum == 0.0) continue;  // NaN: not skipped
      double g = rsum / 2.0;
      double s = 1.0;
      double tot = csum + rsum;
      while (lt_g(csum, g)) {  // dcmpg ifge exit
        s *= 2.0;
        csum *= radix2;
      }
      g = rsum * 2.0;
      while (gt_l(csum, g)) {  // dcmpl ifle exit: NaN EXITS the shrink loop
        s /= 2.0;
        csum /= radix2;
      }
      if (lt_g((csum + rsum) / s, 0.95 * tot)) {
        again = true;
        double f = 1.0 / s;
        for (int k = 0; k < n; ++k) m.ref(i, k) *= f;  // row i, first
        for (int k = 0; k < n; ++k) m.a[idx(i) * idx(n) + idx(k)] *= s;  // column i
      }
    }
  }
}

bool is_hessenberg(const CMat& m) {  // FastMatrix.isHessenberg
  int n = m.n;
  for (int c = 0; c + 2 < n; ++c)
    for (int r = c + 2; r < n; ++r)
      if (m.at(r, c) != 0.0) return false;
  return true;
}

// EigenRoutines.hessenberg (ported for completeness; the JVM never reaches
// it for companion matrices, which are already Hessenberg).
inline void hessenberg_unused(CMat& mat) {
  int n = mat.n;
  for (int mm = 1; mm <= n - 2; ++mm) {
    double scale = 0.0;
    int base = mm * (mm - 1);
    for (int i = mm; i < n; ++i) scale += std::fabs(mat.a[idx(base) + idx(i)]);
    if (scale == 0.0) continue;
    double h = 0.0;
    vd x(static_cast<std::size_t>(n), 0.0);
    std::size_t pos = idx(base) + idx(n - 1);
    for (int i = n - 1; i >= mm; --i) {
      x[idx(i)] = mat.a[pos] / scale;
      pos -= idx(n);
      h += x[idx(i)] * x[idx(i)];
    }
    double g = std::sqrt(h);
    double xp = x[idx(mm)];
    if (xp != 0.0) g = -sign_like(g, xp);
    double f = x[idx(mm)] - g;
    x[idx(mm)] = f;
    mat.a[idx(base) + idx(mm)] = scale * f;
    h -= xp * xp;
    double hh = 0.5 * h;
    if (hh == 0.0) continue;
    vd y(static_cast<std::size_t>(n), 0.0);
    for (int j = 0; j < n; ++j) {
      double g2 = 0.0;
      for (int i = mm; i <= j; ++i) {
        g2 += mat.at(i, j) * x[idx(i)];
        y[idx(i)] = x[idx(i)];
      }
      for (int i = j + 1; i < n; ++i) g2 += mat.at(mm, j) * y[idx(i)];
      y[idx(j)] = g2 / hh;
    }
    f = 0.0;
    for (int j = mm; j < n; ++j) {
      double g2 = 0.0;
      for (int i = mm; i <= j; ++i) {
        g2 += y[idx(i)] * mat.at(i, j);
        x[idx(i)] = y[idx(i)];
      }
      for (int i = j + 1; i < n; ++i) g2 += y[idx(mm)] * mat.at(mm, j);
      x[idx(j)] = g2 / hh;
      f += x[idx(j)] * y[idx(j)];
    }
    double b = f / (hh + hh);
    for (int j = mm; j < n; ++j)
      for (int i = mm; i < n; ++i)
        mat.ref(i, j) -= x[idx(i)] * y[idx(j)] + y[idx(i)] * (x[idx(j)] - b * y[idx(j)]);
  }
  (void)is_hessenberg;
}

// EigenRoutines.hessenbergQR (column-major a, diagonal stride lda = n + 1;
// Support.sign; GeneralEigenSystem sets m_maxiter = 120). Returns the
// eigenvalues; throws jvm_decompose_abort (JVM MatrixException "eig_failed")
// on the maxIter path (never observed in the corpus).
inline std::vector<cplx> hessenberg_qr(CMat a, int maxIter) {
  int n = a.n;
  double nrm = 0.0;
  for (int col = 0; col < n; ++col) {
    int hirow = col < n - 1 ? col + 1 : n - 1;
    for (int pos = col * n; pos <= col * n + hirow; ++pos) nrm += std::fabs(a.a[idx(pos)]);
  }
  std::vector<cplx> ev(static_cast<std::size_t>(n), {0.0, 0.0});
  int k = n - 1;
  int jpos = n * n - 1;
  const int lda = n + 1;
  double accum = 0.0;
  auto aget = [&](int pos) -> double& { return a.a[idx(pos)]; };

  while (k >= 0) {
    int iter = 1;  // slot28, set at the outer head only
    for (;;) {     // inner: the JVM's "continue" (1404) re-runs WITHOUT
      int j = jpos;
      double s = std::fabs(aget(jpos));
      int l = k;
      while (l > 0) {  // deflation search
        j -= lda;
        double hh2 = std::fabs(aget(j));
        double g = hh2 + s;
        if (g == 0.0) g = nrm;  // dcmpl ifne: NaN keeps g
        if (std::fabs(aget(j + 1)) + g != g) {
          s = hh2;
          --l;
        } else {
          aget(j + l) = 0.0;  // verbatim quirky index (harmless; pinned)
          break;
        }
      }
      double q = aget(jpos);
      bool to_outer;
      if (l == k) {  // 1x1 block
        ev[idx(k)] = {q + accum, 0.0};
        --k;
        jpos -= lda;
        to_outer = !(l + 1 < k);
      } else {
        double p = aget(jpos - lda);               // (k-1, k-1)
        double r = aget(jpos - 1) * aget(jpos - n);  // (k-1,k) * (k,k-1)
        if (l == k - 1) {  // 2x2 block
          double x = 0.5 * (p - q);
          double aa = x * x + r;
          double rad = std::sqrt(std::fabs(aa));
          q = q + accum;
          if (!(aa < 0.0)) {  // dcmpl: NaN takes the real branch
            rad = x + sign_like(rad, x);
            ev[idx(k)] = {q + rad, 0.0};
            if (rad != 0.0)
              ev[idx(k - 1)] = {q - r / rad, 0.0};
            else
              ev[idx(k - 1)] = ev[idx(k)];
          } else {
            ev[idx(k)] = {q + x, rad};
            ev[idx(k - 1)] = c_conj(ev[idx(k)]);
          }
          k -= 2;
          jpos -= 2 * lda;
          to_outer = !(l + 1 < k);
        } else {  // general QR iteration
          if (iter == maxIter)
            throw detail::jvm_decompose_abort{"eig_failed"};
          if (iter % 30 == 0) {  // exceptional shift
            accum = accum + q;
            for (int pos = 0, t = 0; t <= k; pos += lda, ++t) aget(pos) -= q;
            double s2 = std::fabs(aget(jpos - n)) + std::fabs(aget(jpos - n - lda));
            p = 0.75 * s2;
            q = p;
            r = (-0.4375 * s2) * s2;
          }
          ++iter;
          int m = k - 2;
          int pos = m * lda;
          int upos = pos + lda;
          int dpos = pos - lda;
          // Slots 20/22/18 (u/v/w) survive the look-ahead loop: the sweep's
          // i == m step SKIPS its load/scale block (@925 if_icmpeq 1015) and
          // reflects with these carry-over values (the (H-s1)(H-s2) first
          // column, already divided by the look-ahead sc). That is also why
          // no a[(m-1)*lda+1] read happens at m == 0.
          double cu = 0.0, cv = 0.0, cw = 0.0;
          while (m >= l) {  // look-ahead
            double f = aget(pos);
            double w0 = q - f;
            double v0 = p - f;
            double u0 = (w0 * v0 - r) / aget(pos + 1) + aget(upos - 1);
            double v1 = aget(upos) - f - w0 - v0;
            double w1 = aget(upos + 1);
            double sc = std::fabs(u0) + std::fabs(v1) + std::fabs(w1);
            u0 /= sc;
            v1 /= sc;
            w1 /= sc;
            cu = u0;
            cv = v1;
            cw = w1;
            if (m == l) break;
            double pp = std::fabs(aget(dpos + 1)) * (std::fabs(v1) + std::fabs(w1));
            double qq = std::fabs(u0) * ((std::fabs(aget(dpos)) + std::fabs(f)) + std::fabs(aget(upos)));
            if (pp + qq == qq) break;
            --m;
            pos -= lda;
            dpos -= lda;
            upos -= lda;
          }
          {  // zero-out bulges
            int i = m;
            int pos2 = m * lda;
            while (i < k - 1) {
              aget(pos2 + 2) = 0.0;      // (i+2, i)
              if (i != m) aget(pos2 - lda + 3) = 0.0;  // (i+2, i-1)
              ++i;
              pos2 += lda;
            }
          }
          {  // Householder sweep (pos starts one diagonal up)
            int i = m;
            int pos2 = (m - 1) * lda;
            while (i < k) {
              double u = cu, v = cv, w = cw;  // slots 20/22/18
              double sc = 0.0;                // slot8 (i == m never reloads)
              if (i != m) {  // @925: the i == m step jumps straight to 1015,
                             // reflecting with the look-ahead carry-over.
                u = aget(pos2 + 1);
                v = aget(pos2 + 2);
                w = 0.0;
                if (i + 1 != k) w = aget(pos2 + 3);
                sc = std::fabs(u) + std::fabs(v) + std::fabs(w);
                if (sc != 0.0) {
                  u /= sc;
                  v /= sc;
                  w /= sc;
                }
              }
              double tau = std::sqrt((u * u + v * v) + w * w);
              double alpha = sign_like(tau, u);
              if (alpha != 0.0) {  // dcmpl ifeq (NaN proceeds)
                if (i == m) {
                  if (m != l) aget(pos2 + 1) = -aget(pos2 + 1);
                } else {
                  aget(pos2 + 1) = -(alpha) * sc;
                }
                u = u + alpha;
                double hv = u / alpha;
                double hu = v / alpha;
                double hw = w / alpha;
                double v2 = v / u;
                double w2 = w / u;
                int pos3 = i * n + i;
                for (int jj = i; jj < k + 1; ++jj, pos3 += n) {  // row update
                  // @1154: s2 = a[pos3] + slot22*a[pos3+1] (+ slot18*a[pos3+2])
                  // where slots 22/18 hold v2 = v/u_new, w2 = w/u_new AFTER
                  // the reassignment at @1118/@1125 (the earlier decode read
                  // them as hv/hw -- corrected 2026-11-21 against the raw
                  // dump; with hv/hw the similarity destroys the trace).
                  double s2 = aget(pos3) + v2 * aget(pos3 + 1);
                  if (i + 1 != k) {
                    s2 += w2 * aget(pos3 + 2);
                    aget(pos3 + 2) = aget(pos3 + 2) - s2 * hw;
                  }
                  aget(pos3 + 1) = aget(pos3 + 1) - s2 * hu;
                  aget(pos3) = aget(pos3) - s2 * hv;
                }
                // @1245: if (k >= i+3) jm = i+3  ELSE jm = k  ==  min(k, i+3)
                int jm = k >= i + 3 ? i + 3 : k;  // column update
                int pr0 = i * n;                  // (verbatim bytecode
                int pr1 = pr0 + n;                //  index arithmetic)
                int pr2 = pr0 + 2 * n;
                for (int j2 = l; j2 < jm + 1; ++j2) {
                  double s2 = hv * aget(pr0 + j2) + hu * aget(pr1 + j2);
                  if (i != k - 1) {
                    s2 += hw * aget(pr2 + j2);
                    aget(pr2 + j2) -= s2 * w2;
                  }
                  aget(pr1 + j2) -= s2 * v2;
                  aget(pr0 + j2) -= s2;
                }
              }
              ++i;
              pos2 += lda;
            }
          }
          to_outer = !(l + 1 < k);
        }
      }
      if (!to_outer) continue;  // inner: same k, iter preserved
      break;                    // outer head: iter resets to 1
    }
  }
  return ev;
}

// GeneralEigenSystem.calc() (pinned): array copy -> balance -> hessenberg
// (only when !isHessenberg(m_std)) -> hessenbergQR, maxIter 120. Shared by
// EV2's roots() and EV1's companion eigen-solve.
inline std::vector<cplx> eigen_pipeline(CMat a) {
  bal(a);
  if (!is_hessenberg(a)) hessenberg_unused(a);
  return hessenberg_qr(a, 120);
}

// EV2 roots(double[]): companion matrix + GeneralEigenSystem (pinned 336-412).
inline std::vector<cplx> ev2_roots(const vd& wts) {
  int m = static_cast<int>(wts.size()) - 1;
  if (m == 0) return {};
  if (m == 1) return {{-(wts[0] / wts[1]), 0.0}};
  CMat a(m);
  double d = wts[idx(m)];
  for (int i = 0; i < m; ++i) a.ref(i, m - 1) = -wts[idx(i)] / d;  // column(m-1)
  for (int t = 0; t + 1 < m; ++t) a.ref(t + 1, t) += 1.0;          // subDiagonal(-1).add(1.0)
  a.ref(0, 1) += 1.0;                                              // add(0, 1, 1.0)
  for (int t = 0; t + 1 < m; ++t) a.ref(t, t + 1) += 1.0;          // subDiagonal(1).add(1.0)
  return eigen_pipeline(a);
}

// ---------------------------------------------------------------------------
// The QR least-squares stack (Reflector.larfg / Householder2 /
// QRDecomposition.applyQt / UpperTriangularMatrix / QRLeastSquaresSolver /
// LeastSquaresDivision), column-major ld = m, pinned byte-exact.
// ---------------------------------------------------------------------------

// DataPointer.norm2(n): scaled accumulation x*sqrt(1+sum(r^2)).
inline double dp_norm2(const double* p, int n, int inc) {
  if (n < 1) return 0.0;
  if (n == 1) return std::fabs(p[0]);
  double x = 0.0, ss = 1.0;
  for (int t = 0; t < n; ++t) {
    double xi = p[idx(t) * idx(inc)];
    if (xi != 0.0) {
      double ax = std::fabs(xi);
      if (ax > x) {
        double r = x / ax;
        ss = 1.0 + ss * r * r;
        x = ax;
      } else {
        double r = ax / x;
        ss = ss + r * r;
      }
    }
  }
  return x * std::sqrt(ss);
}

// Reflector.larfg over the window (storage + xstart, length n).
struct Larfg {
  double x0;
  double beta;
  double alpha;
};

inline Larfg larfg(std::vector<double>& storage, int xstart, int n) {
  Larfg o{0.0, 0.0, 0.0};
  if (n == 0) return o;  // JVM: void larfg() early return (beta stays 0)
  double x0 = storage[idx(xstart)];  // x0() reads the head slot
  double* v = storage.data() + xstart + 1;  // v pointer (tail window)
  int nm1 = n - 1;
  auto x_norm2 = [&]() { return dp_norm2(storage.data() + xstart, n, 1); };
  double nrm = x_norm2();  // x.norm2(n) INCLUDES the head slot
  constexpr double eps = kFeEps;                    // Constants.getEpsilon()
  double ssfmin = 2.2250738585072014e-308 / eps;    // getSafeMin()/epsilon
  int it = 0;
  if (nrm != 0.0 && lt_g(nrm, ssfmin)) {  // dcmpl ifeq / dcmpg ifge
    double sc = 1.0 / ssfmin;             // computed ONCE (bytecode @74)
    for (;;) {
      for (int i = 0; i < nm1; ++i) v[i] *= sc;  // v.mul(n-1, sc)
      x0 *= sc;
      nrm *= sc;
      if (!lt_g(nrm, ssfmin)) break;  // @105 dcmpg ifge -> recompute
      ++it;                           // @108: counted only on a REPEAT
      if (it < 4) continue;           // @114 if_icmplt -> scale again
      break;                          // @117 falls through too
    }
    nrm = x_norm2();  // @117: the norm of the SCALED vector is always
                      // recomputed when the scaling block ran (it stays
                      // scaleCount-1; the JVM's undo loop below is pinned
                      // verbatim, off-by-one quirk included)
  }
  if (lt_g(x0, 0.0)) nrm = -nrm;                         // @128 dcmpg
  for (int i = 0; i < it; ++i) nrm *= ssfmin;            // @139
  o.beta = (nrm == 0.0) ? 0.0 : nrm / (nrm + x0);        // @162 dcmpl ifne
  for (int i = 0; i < nm1; ++i) v[i] /= nrm;             // @184 (NO guard:
                                                         // div-by-zero -> Inf)
  o.alpha = -nrm;                                        // @202
  // @193-199 writes x0 = 1 + x0/nrm back into the head slot; Householder2
  // overwrites the slot with alpha, so the value is unobservable here.
  return o;
}

// Reflector.lapply over the trailing window (columns of the shrunk matrix).
void lapply(std::vector<double>& qr, int xstart, int n, double beta, int win_start,
            int win_cols, int ld) {
  int vend = xstart + n;
  for (int j = 0; j < win_cols; ++j) {
    int base = win_start + j * ld;
    double tau = qr[idx(base)] / beta;
    {
      int q = base + 1;
      for (int i = xstart + 1; i < vend; ++i, ++q) tau += qr[idx(q)] * qr[idx(i)];
    }
    if (tau == 0.0) continue;
    qr[idx(base)] -= tau;
    tau = -tau * beta;
    {
      int q = base + 1;
      for (int i = xstart + 1; i < vend; ++i, ++q) qr[idx(q)] += tau * qr[idx(i)];
    }
  }
}

struct QRDecomp {
  std::vector<double> qr;  // column-major, ld = m
  std::vector<double> beta;
  int m = 0, n = 0;

  // QRDecomposition.applyQt(double[]) (ascending k).
  void apply_qt(vd& x) const {
    if (static_cast<int>(x.size()) != m)
      throw detail::jvm_decompose_abort{"m_err_dim"};
    for (int k = 0, pos = 0; k < n; ++k, pos += m + 1) {
      double b = beta[idx(k)];
      if (b == 0.0) continue;
      double tau = x[idx(k)] / b;
      for (int i = k + 1, q = pos + 1; i < m; ++i, ++q) tau += x[idx(i)] * qr[idx(q)];
      if (tau == 0.0) continue;
      x[idx(k)] -= tau;
      tau = -tau * b;
      for (int i = k + 1, q = pos + 1; i < m; ++i, ++q) x[idx(i)] += tau * qr[idx(q)];
    }
  }

  // rawR(): extract(0,n,0,n).deepClone() + toUpper.
  vd raw_r() const {
    vd r = qr;  // full storage clone (ld = m preserved)
    for (int col = 1; col < n; ++col)
      for (int row = col + 1; row < n; ++row)
        r[idx(row) + idx(col) * idx(m)] = 0.0;
    return r;
  }
};

// Householder2.decompose.
inline QRDecomp householder2(std::vector<double> a, int m, int n) {
  QRDecomp q;
  q.qr = std::move(a);
  q.m = m;
  q.n = n;
  q.beta.assign(static_cast<std::size_t>(n), 0.0);
  int kmax = std::min(m, n);
  int pos = 0;
  for (int k = 0; k < kmax; ++k) {
    Larfg lr = larfg(q.qr, pos, m - k);
    int cols = n - k - 1;
    if (cols > 0) {
      // window.bhshrink() always advances; lapply runs only when beta != 0
      // (the pinned quirk).
      if (lr.beta != 0.0) lapply(q.qr, pos, m - k, lr.beta, pos + m, cols, m);
      q.qr[idx(pos)] = lr.alpha;
    } else {
      q.qr[idx(pos)] = lr.alpha;
    }
    q.beta[idx(k)] = lr.beta;
    pos += m + 1;
  }
  return q;
}

struct LSResult {
  bool ok = false;  // false = MatrixException path
  vd coeff;
  double ssq_err = 0.0;
};

// QRLeastSquaresSolver.leastSquares (tol = 0) on an existing decomposition.
inline LSResult least_squares(const QRDecomp& qr, const vd& y) {
  LSResult r;
  vd R = qr.raw_r();
  int rank = 0;  // fastRank(rawR, 0): |d| > 0 (NaN not counted)
  for (int i = 0; i < qr.n; ++i) {
    double d = R[idx(i) + idx(i) * idx(qr.m)];
    if (std::fabs(d) > 0.0) ++rank;
  }
  vd b = y;
  try {
    qr.apply_qt(b);
  } catch (const detail::jvm_decompose_abort&) {
    r.ok = false;
    return r;
  }
  // solveUx(rawR[0..rank,0..rank], b, tol = 0): backward substitution.
  for (int i = rank - 1; i >= 0; --i) {
    double d = R[idx(i) + idx(i) * idx(qr.m)];
    if (!(std::fabs(d) > 0.0)) {
      // |xi| < 0 never holds; the dcmpl tol=0 comparison jumps (x = 0) only
      // when |xi| is NaN. Everything else throws m_err_sing.
      if (!is_nan(b[idx(i)])) {
        r.ok = false;  // MatrixException("m_err_sing")
        return r;
      }
      b[idx(i)] = 0.0;
      continue;
    }
    double v = b[idx(i)] / d;
    b[idx(i)] = v;
    for (int j = i - 1, up = static_cast<int>(idx(i) + idx(i) * idx(qr.m)) - 1; j >= 0;
         --j, --up)
      b[idx(j)] -= v * R[idx(up)];
  }
  double ssq = 0.0;
  for (int i = rank; i < qr.m; ++i) ssq += b[idx(i)] * b[idx(i)];
  if (rank != qr.n) {
    r.coeff.assign(static_cast<std::size_t>(qr.n), 0.0);
    for (int i = 0; i < rank; ++i) r.coeff[idx(i)] = b[idx(i)];
  } else {
    r.coeff.assign(b.begin(), b.begin() + rank);
  }
  r.ssq_err = ssq;
  r.ok = true;
  return r;
}

// LeastSquaresDivision.divide(p, q); err = ssqErr / q.length (pinned). The
// JVM keeps the quotient field even when solveSquares failed; getQuotient()
// then returns whatever was set (null on the never-divided path -> NPE), so
// the mirror keeps the same split: false keeps `quotient' untouched.
struct LsDiv {
  Poly quotient;
  bool have_quotient = false;
  double err = 0.0;
};

inline LsDiv ls_divide(const Poly& p, const Poly& q) {
  LsDiv out;
  int m = static_cast<int>(p.c.size());
  int n = static_cast<int>(q.c.size());
  if (n > m) return out;  // divide() returns false, quotient stays null
  int ncol = m - n + 1;
  std::vector<double> A(static_cast<std::size_t>(m) * static_cast<std::size_t>(ncol), 0.0);
  for (int j = 0; j < ncol; ++j)
    for (int i = 0; i < n; ++i) A[idx(j + i) + idx(j) * idx(m)] = q.c[idx(i)];
  QRDecomp qr = householder2(std::move(A), m, ncol);
  LSResult r = least_squares(qr, p.c);
  if (!r.ok) return out;  // MatrixException -> divide() false, quotient set
  out.quotient = p_of(r.coeff);  // getQuotient(): Polynomial.of(coefficients)
  out.have_quotient = true;
  out.err = r.ssq_err / static_cast<double>(n);
  return out;
}

// LeastSquaresDivision.isExact: err < sqrt(EPSILON).
bool ls_exact(double err) { return err < std::sqrt(kFeEps); }

// ---------------------------------------------------------------------------
// Polynomial$SimplifyingTool + BackFilter$SimplifyingTool (the pm /
// simplifyAr common-factor extraction; the LS-exact branch is checked
// first, then the numeric roots branch - both transcribed verbatim,
// including the shared-tail double divide the bytecode performs).
// ---------------------------------------------------------------------------
struct PolySimp {
  bool ret = false;
  Poly common = p_of(vd{1.0});
  bool has_lr = false;
  Poly left, right;
};

// Polynomial.evaluateAt(Complex) (ComplexType Horner).
inline cplx eval_complex(const Poly& p, cplx z) {
  cplx acc{0.0, 0.0};
  for (int i = p.degree(); i >= 0; --i) acc = c_add(c_times(acc, z), c_real(p.c[idx(i)]));
  return acc;
}

inline PolySimp poly_simplify_private(const Poly& big, const Poly& small) {
  PolySimp t;
  if (small.c.size() == 1) return t;  // false, common = ONE
  LsDiv dv = ls_divide(big, small);
  if (dv.have_quotient && ls_exact(dv.err)) {  // simplifyExact
    t.has_lr = true;
    t.left = dv.quotient;
    t.right = p_of(vd{1.0});
    t.common = small;
    t.ret = true;
    return t;
  }
  std::vector<cplx> zroots;
  try {
    // 24-C.1 cluster 3 RE-ADJUDICATION (was SECTION 15 divergence (a)/(b)):
    // the pinned bytecode takes roots of ARG2 (the smaller polynomial) and a
    // NaN imaginary part REFRESHES f2 (dcmpl ifle on im <= 0 falls through
    // with the stale f2; NaN gets dcmpl -> 1 and refreshes it). Both pinned
    // here; the decompose corpus never enters this path (its 49/1669 replay
    // re-verified byte-identical after the fix).
    zroots = detail::p_roots(small, "simplify");  // Polynomial.roots() on ARG2
  } catch (const std::invalid_argument&) {
    throw;  // envelope violation, never swallowed
  }
  vd f1 = {0.0, 1.0};   // Coefficients.fromDegree(1)
  vd f2 = {0.0, 0.0, 1.0};  // Coefficients.fromDegree(2)
  Poly common = p_of(vd{1.0});
  for (const cplx& z : zroots) {
    if (lt_g(c_abs(eval_complex(big, z)), 1.0E-12)) {  // dcmpg: NaN skips
      double re = z.re, im = z.im;
      if (im == 0.0) {
        f1[0] = -re;
        Poly q1 = p_raw(pl_divide(big.c, f1));
        common = p_raw(pl_times(common.c, q1.c));  // new Polynomial(raw)
      }
      // 149: goto 184 - the shared f2 tail ALSO runs after the real branch;
      // im <= 0.0 (dcmpl ifle) keeps the stale f2, NaN im refreshes it.
      if (!(im <= 0.0)) {
        f2[0] = re * re + im * im;
        f2[1] = -2.0 * re;
      }
      Poly q2 = p_raw(pl_divide(big.c, f2));
      common = p_raw(pl_times(common.c, q2.c));
    }
  }
  if (common.degree() > 0) {
    t.has_lr = true;
    t.left = big;
    t.right = p_quotient(small, common);  // small.divide(common)
    t.common = common;
    t.ret = true;
  } else {
    t.common = common;
  }
  return t;
}

PolySimp poly_simplify(const Poly& p, const Poly& q) {  // public simplify
  PolySimp t;
  if (p.c.size() >= q.c.size()) {
    t = poly_simplify_private(p, q);
  } else {
    t = poly_simplify_private(q, p);
    std::swap(t.left, t.right);  // the JVM swaps the (possibly unset) fields
  }
  t.ret = t.common.degree() > 0;
  return t;
}

// ---------------------------------------------------------------------------
// BackFilter (one-sided filter = the Polynomial itself) and its
// SimplifyingTool.
// ---------------------------------------------------------------------------
struct BFilter {
  Poly p;
  int length() const { return p.degree() + 1; }
  int degree() const { return p.degree(); }
  static BFilter one() { return BFilter{p_of(vd{1.0})}; }
};

BFilter bf_times(const BFilter& x, const BFilter& y) {  // BackFilter.times
  return BFilter{p_times(x.p, y.p)};
}

struct BackSimp {
  bool ret = false;
  BFilter common, left, right;
};

BackSimp back_simplify(const BFilter& a, const BFilter& b) {  // tool
  BackSimp out;
  if (a.degree() == 0 || b.degree() == 0) return out;
  double a0 = a.p.c[0], b0 = b.p.c[0];
  PolySimp t = poly_simplify(a.p, b.p);
  if (!t.ret) return out;
  Poly p = t.left, q = t.right;
  p = p_times_d(p, a0 / p.c[0]);
  q = p_times_d(q, b0 / q.c[0]);
  Poly common = p_times_d(t.common, 1.0 / t.common.c[0]);  // divide(get(0))
  out.ret = true;
  out.common = BFilter{std::move(common)};
  out.left = BFilter{std::move(p)};
  out.right = BFilter{std::move(q)};
  return out;
}

// ---------------------------------------------------------------------------
// Shared removeUnitRoots of the two eigenvalue decomposers (the JVM source
// is literally identical in EigenValuesDecomposer and EigenValuesDecomposer2:
// UnitRootsSolver factorize -> UnitRoots.sqrt -> bf = BackFilter(h) -> LS
// divide by h^2 (return discarded) -> getQuotient (unconditional read: a
// failed divide NPEs exactly like the JVM) -> symmetrize -> Polynomial.of;
// bf defaults to BackFilter.ONE).
// ---------------------------------------------------------------------------
inline Poly ev_remove_unit_roots(const Poly& p_in, std::optional<BFilter>& bf) {
  UnitRoots ur = ur_solver_factorize(p_in);
  if (ur.roots_count() == 0) {  // solver.factorize(): getRootsCount() != 0
    if (!bf) bf = BFilter::one();
    return p_in;
  }
  Poly p = p_in;
  std::optional<UnitRoots> h = ur.sqrt_();
  if (h.has_value()) {
    Poly hp = h->as_polynomial();
    Poly hh = p_times(hp, hp);
    bf = BFilter{hp};
    LsDiv dv = ls_divide(p, hh);  // divide()'s return value is discarded
    if (!dv.have_quotient) throw JvmNpe{};  // getQuotient() null -> NPE
    vd a = dv.quotient.c;
    int d = dv.quotient.degree();
    int n = d / 2;
    for (int i = 0; i < n; ++i) {
      double m = (a[idx(i)] + a[idx(d - i)]) / 2.0;
      a[idx(i)] = m;
      a[idx(d - i)] = m;
    }
    p = p_of(std::move(a));  // Polynomial.of (trims)
  }
  if (!bf) bf = BFilter::one();
  return p;
}

// ---------------------------------------------------------------------------
// EigenValuesDecomposer2 (the SymmetricFilter factorizer installed first).
// decompose() returns false exactly where the JVM does: any JVM exception
// (NPE, AIOOBE, MatrixException from hessenbergQR, ...) is caught -> false.
// A DECLINE is not final: the chain continues with EV1 (and further), so
// factorize_sf must not stop here.
// ---------------------------------------------------------------------------
struct Ev2 {
  std::optional<BFilter> bf;  // null == the JVM null field
  double fac = 0.0;

  Poly remove_unit_roots(const Poly& p_in) {
    return ev_remove_unit_roots(p_in, bf);
  }

  bool decompose(const SymF& sf) {
    try {
      bf = std::nullopt;
      fac = 0.0;
      double w0 = sf_weight(sf, 0);
      if (w0 < 0.0) return false;  // dcmpg
      if (sf.a.size() == 1) {
        bf = BFilter::one();
      } else {
        Poly p = p_of(sf_weights_array(sf));  // Polynomial.of(palindrome)
        p = remove_unit_roots(p);
        int n = p.degree() / 2;
        if (n > 0) {
          vd wts(p.c.begin() + n, p.c.begin() + n + n + 1);  // extract(n, n+1)
                                                              // == c[n..2n]
          std::vector<cplx> zz = ev2_roots(wts);
          std::vector<cplx> ev(zz.size(), {0.0, 0.0});
          vd ur(zz.size() / 2, 0.0);
          std::size_t cnt = 0, eidx = 0;
          for (const cplx& z : zz) {
            if (z.im == 0.0) {  // dcmpl ifne
              double re = z.re;
              double dsc = re * re - 4.0;
              if (lt_g(dsc, -1.0E-9)) {  // dcmpg
                cplx qc{re / 2.0, std::sqrt(-dsc) / 2.0};
                std::size_t j = 0;
                while (j < cnt) {
                  if (lt_g(std::fabs(re - ur[idx(j)]), 1.0E-6)) {  // dcmpg
                    ur[idx(j)] = 0.0;
                    break;
                  }
                  ++j;
                }
                if (j == cnt) {
                  if (cnt >= ur.size()) throw Aioobe{};  // new double[len/2]
                  ur[idx(cnt++)] = re;
                  if (eidx + 1 >= ev.size()) throw Aioobe{};  // JVM array
                                                              // bounds
                  ev[idx(eidx++)] = qc;
                  ev[idx(eidx++)] = c_conj(qc);
                }
              } else {
                if (lt_g(dsc, 0.0)) {  // dcmpg: NaN skips BOTH clamps
                  dsc = 0.0;
                  re = (re < 0.0) ? -2.0 : 2.0;  // dcmpg: NaN -> +2
                }
                double sq = std::sqrt(dsc);
                if (eidx >= ev.size()) throw Aioobe{};
                if (re < 0.0)  // dcmpg: NaN -> the + branch
                  ev[idx(eidx++)] = c_cart((re - sq) / 2.0, 0.0);
                else
                  ev[idx(eidx++)] = c_cart((re + sq) / 2.0, 0.0);
              }
            } else if (gt_l(z.im, 0.0)) {  // dcmpl ifle: im < 0 / NaN drop
              cplx s = c_sqrt(c_minus(c_times(z, z), c_real(4.0)));
              cplx cp{(z.re + s.re) / 2.0, (z.im + s.im) / 2.0};
              cplx cm{(z.re - s.re) / 2.0, (z.im - s.im) / 2.0};
              cplx sel = gt_l(c_abs(cp), c_abs(cm)) ? cp : cm;  // dcmpl
              if (eidx + 1 >= ev.size()) throw Aioobe{};
              ev[idx(eidx++)] = sel;
              ev[idx(eidx++)] = c_conj(sel);
            }
          }
          if (eidx != ev.size()) throw JvmNpe{};  // lejaOrder hits a null
          std::vector<cplx> full(ev.size());
          for (std::size_t i = 0; i < ev.size(); ++i) full[idx(i)] = ev[idx(i)];
          leja_order(full);
          Poly pp = p_from_roots(full);
          pp = p_times_d(pp, 1.0 / pp.c[0]);
          bf = bf_times(*bf, BFilter{pp});
        }
      }
      double ssq = 0.0;
      for (double x : bf->p.c) ssq += x * x;
      fac = w0 / ssq;
      return true;
    } catch (const std::invalid_argument&) {
      throw;  // envelope violation, never swallowed
    } catch (...) {
      return false;  // JVM: catch (Exception e) { return false; }
    }
  }
};

// ---------------------------------------------------------------------------
// internal.toolkit.base.core.math.linearfilters.EigenValuesDecomposer (EV1):
// the factorizer tried after EV2 declines. Bytecode dump
// jdprobe/m524/EigenValuesDecomposer.txt (decompose 27-208,
// removeUnitRoots 210-304). Unlike EV2 it keeps the FULL-degree companion
// (subDiagonal(-1).drop(0,1) set to 1, last column -a[i]/a[n] with a[n] the
// top coefficient, and (0, n) = 1.0), then filters the eigenvalues to
// |z| > 1 (dropping |z| <= 1e-9 first), leja-orders them and builds the B
// filter from the roots. corpus: supplies the 5 components where EV2
// overflows its own ev[] array (v_simplify r c1, d_d2_decomp r c1 - the
// rest are refused by the validator before any factorizer runs).
// ---------------------------------------------------------------------------
struct Ev1 {
  std::optional<BFilter> bf;
  double fac = 0.0;

  bool decompose(const SymF& sf) {
    try {
      bf = std::nullopt;
      fac = 0.0;
      double w0 = sf_weight(sf, 0);
      if (w0 < 0.0) return false;  // dcmpg ifge
      if (sf.a.size() == 1) {
        bf = BFilter::one();
      } else {
        vd a = sf_weights_array(sf);
        for (double& x : a) x /= w0;  // @53-73: in-place, EVERY entry
        Poly p = p_of(std::move(a));
        p = ev_remove_unit_roots(p, bf);
        int n = p.degree();  // FULL degree (EV2 used degree/2)
        if (n > 0) {
          vd pa = p.c;
          CMat mat(n + 1);
          for (int t = 0; t + 1 < n; ++t)
            mat.ref(t + 1, t) = 1.0;  // subDiagonal(-1).drop(0, 1).set(1)
          // column(n-1).drop(0,1) = rows 0..n-1 (drop(0 front, 1 back)) <-
          // setAY(-1/a[n], DataBlock.of(a, 0, n)).
          double top = -1.0 / pa[idx(n)];
          for (int i = 0; i < n; ++i) mat.ref(i, n - 1) = top * pa[idx(i)];
          mat.ref(0, n) = 1.0;  // FastMatrix.set(0, n, 1)
          std::vector<cplx> roots = eigen_pipeline(mat);
          std::vector<cplx> ev(roots.size() / 2, {0.0, 0.0});
          std::size_t eidx = 0;
          for (const cplx& z : roots) {
            if (!gt_g(c_abs(z), 1.0E-9)) continue;  // equals(ZERO, 1e-9)
            if (!gt_l(c_abs(z), 1.0)) continue;     // |z| > 1 (dcmpl)
            if (eidx >= ev.size()) throw Aioobe{};  // JVM array bounds
            ev[idx(eidx++)] = z;
          }
          if (eidx != ev.size()) throw JvmNpe{};  // nulls reach lejaOrder
          leja_order(ev);
          Poly pp = p_from_roots(ev);
          pp = p_times_d(pp, 1.0 / pp.c[0]);
          bf = bf_times(*bf, BFilter{pp});
        }
      }
      double ssq = 0.0;
      for (double x : bf->p.c) ssq += x * x;
      fac = w0 / ssq;
      return true;
    } catch (const std::invalid_argument&) {
      throw;  // envelope violation, never swallowed
    } catch (...) {
      return false;  // JVM: catch (Exception e) { return false; }
    }
  }
};

// SymmetricFilter$Factorization.
struct Factorization {
  BFilter bf;
  double scaling;
};

// SymmetricFilter.factorize(): the degree-0 shortcut first (negative and
// zero constants are legal), then the factorizer chain. Corpus-verified
// over every distinct symmetric filter (150 of them): 124 fall to EV2, 5
// to EV1, 1 to the fast decomposer, 20 to none. The single fast-decomposer
// filter lives in a case the validator refuses before any factorizer runs
// (v_airline r c1), so the mirror chains EV2 -> EV1 and stops; robust is
// never reached.
inline std::optional<Factorization> factorize_sf(const SymF& sf) {
  if (sf.a.size() == 1) return Factorization{BFilter::one(), sf.a[0]};
  Ev2 ev;
  if (ev.decompose(sf)) return Factorization{*ev.bf, ev.fac};
  Ev1 ev1;
  if (ev1.decompose(sf)) return Factorization{*ev1.bf, ev1.fac};
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// jdplus ...core.arima.ArimaModel (the symmetric and one-sided variants;
// the lazily cached sar/derivedsma/factorization slots are recomputed per
// call - they are pure functions of the stored state, and the JVM never
// caches a null factorization).
// ---------------------------------------------------------------------------
struct AModel {
  Poly ar, delta;
  bool symmetric = true;
  SymF sma;       // set when symmetric
  vd ma1;         // set when !symmetric (the one-sided MA)
  double var = 0.0;
  bool has_ma = false;  // the JVM `ma' field is non-null

  const Poly& stationary_ar() const { return ar; }
  const Poly& non_stationary_ar() const { return delta; }

  BFilter get_ma() const;               // throws MaError on null factorization
  double get_iv() const;                // throws JvmNpe on null factorization
  SymF symmetric_ma() const {
    if (symmetric) return sma;
    return sym_convolution(ma1, var);
  }
  SymF symmetric_ar() const {
    return sym_convolution(p_times(ar, delta).c, 1.0);
  }
  bool is_white_noise() const {
    if (ar.c.size() > 1) return false;
    if (delta.c.size() > 1) return false;
    if (has_ma && ma1.size() > 1) return false;
    if (symmetric && sma.a.size() > 1) return false;
    return gt_l(get_iv(), 0.0);  // dcmpl ifle: NaN -> false
  }
  bool is_null() const {
    if (delta.c.size() != 1) return false;
    if (symmetric) return sf_is_null(sma);
    return std::fabs(var) < 1.0E-6;  // dcmpg: NaN -> false
  }
  detail::Spectrum spectrum() const {
    detail::Spectrum s;
    s.num = symmetric_ma();
    s.den = symmetric_ar();
    return s;
  }
};

// ArimaModel.NULL = whiteNoise(0.0) (static initializer). The isNull()
// formula reproduces the singleton check for every model we build.
inline const AModel& am_null() {
  static const AModel n = [] {
    AModel a;
    a.ar = p_of(vd{1.0});
    a.delta = p_of(vd{1.0});
    a.symmetric = false;
    a.has_ma = true;
    a.ma1 = {1.0};
    a.var = 0.0;
    return a;
  }();
  return n;
}

inline AModel am_white_noise(double d) {
  AModel a = am_null();
  a.var = d;
  return a;
}

inline BFilter AModel::get_ma() const {
  if (has_ma) return BFilter{p_raw(vd(ma1))};
  std::optional<Factorization> fac = factorize_sf(sma);
  if (!fac) throw MaError{};  // ArimaException("invalid ARIMA")
  return fac->bf;
}

inline double AModel::get_iv() const {
  if (has_ma) return var;
  std::optional<Factorization> fac = factorize_sf(sma);
  if (!fac) throw JvmNpe{};  // reading scaling of a null Factorization
  return fac->scaling;
}

// ArimaModel.minus(double) (pinned).
inline AModel am_minus_d(const AModel& a, double d) {
  if (a.is_white_noise()) {
    double diff = a.var - a.get_iv();  // Math.abs(var - getIV()) < 1e-6
    if (lt_g(std::fabs(diff), 1.0E-6)) return am_null();
  }
  SymF sma = a.symmetric_ma();
  SymF sar = a.symmetric_ar();
  SymF newsma = sf_minus(sma, sf_multiply_d(d, sar));
  AModel r;
  r.ar = a.ar;
  r.delta = a.delta;
  r.symmetric = true;
  r.sma = std::move(newsma);
  return r;
}

// ArimaModel.plus(double) (pinned).
inline AModel am_plus_d(const AModel& a, double d) {
  if (a.is_white_noise()) return am_white_noise(d + a.get_iv());
  SymF sma = a.symmetric_ma();
  SymF sar = a.symmetric_ar();
  SymF newsma = sf_plus(sma, sf_multiply_d(d, sar));
  AModel r;
  r.ar = a.ar;
  r.delta = a.delta;
  r.symmetric = true;
  r.sma = std::move(newsma);
  return r;
}

// ArimaModel.simplifyAr() (pinned): the BackFilter tool on (ar, getMa()).
// The ArimaException from getMa()/the tool is caught -> this; anything else
// (NPE from getIV) propagates exactly like the JVM.
inline AModel am_simplify_ar(const AModel& a) {
  try {
    BFilter ma = a.get_ma();
    BackSimp t = back_simplify(BFilter{a.ar}, ma);
    if (!t.ret) return a;
    AModel r;
    r.ar = t.left.p;
    r.delta = a.delta;
    r.symmetric = false;
    r.has_ma = true;
    r.ma1 = t.right.p.c;
    r.var = a.get_iv();  // NPE not caught
    return r;
  } catch (const MaError&) {
    return a;
  }
}

// ArimaModel.pm(a, b, add, trim) (pinned, transcribed statement by
// statement). Returns NULL (the singleton) when newsma.isNull().
inline AModel am_pm(const AModel& a, const AModel& b, bool add, bool trim) {
  if (b.is_white_noise()) return add ? am_plus_d(a, b.get_iv()) : am_minus_d(a, b.get_iv());
  if (a.is_white_noise() && add) return am_plus_d(b, a.get_iv());  // add(a.getIV(), b)
  BFilter ar1{a.ar}, ar2{b.ar};
  BFilter delta12;
  {
    BackSimp t = back_simplify(ar1, ar2);
    if (t.ret) {
      delta12 = bf_times(ar1, t.right);
      ar1 = t.left;
      ar2 = t.right;
    } else {
      delta12 = bf_times(ar1, ar2);
    }
  }
  BFilter d1{a.delta}, d2{b.delta}, ar;
  {
    BackSimp t = back_simplify(d1, d2);
    if (t.ret) {
      ar = bf_times(d1, t.right);
      ar1 = bf_times(ar1, t.left);
      ar2 = bf_times(ar2, t.right);
    } else {
      ar = bf_times(d1, d2);
      ar1 = bf_times(ar1, d1);
      ar2 = bf_times(ar2, d2);
    }
  }
  SymF conv1 = sym_convolution(ar1.p.c, 1.0);
  SymF conv2 = sym_convolution(ar2.p.c, 1.0);
  SymF sma1 = a.symmetric_ma();
  SymF sma2 = b.symmetric_ma();
  SymF newsma = add ? sf_plus(sf_times(sma1, conv2), sf_times(sma2, conv1))
                    : sf_minus(sf_times(sma1, conv2), sf_times(sma2, conv1));
  if (sf_is_null(newsma)) return am_null();
  AModel r;
  r.ar = delta12.p;
  r.delta = ar.p;
  r.symmetric = true;
  r.sma = std::move(newsma);
  // JVM pm tail (pinned ArimaModel.txt 335-350): `return trim
  // ? model.simplify() : model.simplifyAr();' - simplifyAr() ALWAYS runs, the
  // trim flag only upgrades it to simplify() (simplifyAr + simplifyUr).
  // trim is false on every path reachable from 24-B (both minus() and plus()
  // pass it), so the simplify() branch is documented-unreachable here.
  (void)trim;
  return am_simplify_ar(r);
}

// ArimaModel.plus(b, trim) / minus(b, trim).
inline AModel am_plus(const AModel& a, const AModel& b, bool trim) {
  return am_pm(a, b, true, trim);
}

// ---------------------------------------------------------------------------
// jdplus ...core.ucarima.UcarimaModel.
// ---------------------------------------------------------------------------
struct UcarimaModel {
  // The outer model (the SarimaModel on every corpus path): dumped straight
  // from the raw getters; never an ArimaModel inside the envelope.
  vd m_sar, m_ur, m_ma;
  double m_var = 1.0;
  std::vector<AModel> comps;
};

// UcarimaModel.isValid(): any getIV() < 0 -> false; ArimaException -> false;
// NPE propagates (the probe would print "excep").
inline bool um_is_valid(const UcarimaModel& um) {
  try {
    for (const AModel& c : um.comps) {
      double iv = c.get_iv();
      if (lt_g(iv, 0.0)) return false;  // dcmpg: NaN continues
    }
    return true;
  } catch (const MaError&) {
    return false;
  }
}

// UcarimaModel.getComplement(i) (pinned UcarimaModel.txt 110, 24-C.1 cluster
// 3): NULL when the model has a single component; else the fold
// NULL.plus(c, false) over the components with j != i and !c.isNull(). The
// JVM caches the fold in complements[i] (synchronized) - a pure function of
// the stored components, so the mirror recomputes per call. The JVM reads
// complements[i] BEFORE the fold (bytecode 13-18): out-of-range i throws the
// AIOOBE even when the fold would not touch index i (mirrored).
inline AModel um_get_complement(const UcarimaModel& um, int i) {
  int n = static_cast<int>(um.comps.size());
  if (n <= 1) return am_null();  // if_icmpgt: length > 1 only
  if (i < 0 || i >= n)
    throw std::out_of_range("seats_ucarima_model: getComplement (AIOOBE parity)");
  AModel r = am_null();
  for (int j = 0; j < n; ++j) {
    if (j == i) continue;
    const AModel& c = um.comps[idx(j)];
    if (c.is_null()) continue;
    r = am_plus(r, c, false);
  }
  return r;
}

// UcarimaModel.setVarianceMax(-1, false) (pinned, including the branch
// structure: comps[n].plus(m) runs ONLY when !(m < 0), NaN takes the plus
// branch, and the noisy branch would replace the model - unreachable in the
// envelope, still transcribed for completeness).
inline std::optional<UcarimaModel> um_set_variance_max(const UcarimaModel& um, int n, bool noisy) {
  double m = 0.0;
  if (um.comps.empty()) return um;
  std::vector<AModel> comps;
  if (n < 0) {
    comps = um.comps;
    comps.push_back(am_null());
    n = static_cast<int>(um.comps.size());
  } else {
    comps = um.comps;
  }
  for (int i = 0; i < static_cast<int>(comps.size()); ++i) {
    if (i == n) continue;
    const AModel& c = comps[idx(i)];
    detail::MinResult mr = detail::minimize_spectrum(c.spectrum());
    if (mr.min == 0.0) continue;  // dcmpl ifeq: NaN proceeds
    m += mr.min;
    comps[idx(i)] = am_minus_d(c, mr.min);
  }
  UcarimaModel out = um;
  if (lt_g(m, 0.0)) {  // dcmpg: NaN takes the plus branch
    if (noisy) {
      return std::nullopt;  // the Noisy branch is refused before we get here
    }
    return std::nullopt;  // JVM: return null
  }
  comps[idx(n)] = am_plus_d(comps[idx(n)], m);
  out.comps = std::move(comps);
  return out;
}

// UcarimaModel.compact(p, s) (pinned) with the verify=false Builder.
inline UcarimaModel um_compact(const UcarimaModel& um, int p, int s) {
  UcarimaModel out;
  out.m_sar = um.m_sar;
  out.m_ur = um.m_ur;
  out.m_ma = um.m_ma;
  out.m_var = um.m_var;
  for (int i = 0; i < p; ++i) out.comps.push_back(um.comps[idx(i)]);
  AModel comp = um.comps[idx(p)];
  for (int i = 1; i < s; ++i) comp = am_plus(comp, um.comps[idx(p + i)], false);
  out.comps.push_back(comp);
  for (int i = p + s; i < static_cast<int>(um.comps.size()); ++i)
    out.comps.push_back(um.comps[idx(i)]);
  return out;
}

// SeatsModel.compactUcarimaModel(false, false): compact(2, 2) iff count > 3.
inline UcarimaModel seats_compact_ucarima_model(const UcarimaModel& um) {
  if (static_cast<int>(um.comps.size()) > 3) return um_compact(um, 2, 2);
  return um;
}

}  // namespace seats_model
}  // namespace detail
}  // namespace kalos
