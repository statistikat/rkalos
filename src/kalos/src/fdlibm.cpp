#include "kalos/detail/fdlibm.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

// Operation-for-operation transcription of java.lang.FdLibm.Exp.compute and
// java.lang.FdLibm.Log.compute as pinned in the reference runtime (bytecode
// dumps: FdExp.txt, FdLog.txt). Every double expression below keeps the
// operand order of the corresponding dadd/dsub/dmul/ddiv sequence, including
// the redundant-looking groupings; the JVM's int arithmetic is 32-bit wrapping,
// mirrored here on uint32_t with the signed comparisons the bytecode makes.

namespace kalos::detail {
namespace {

uint32_t hi_word(const double x) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof(bits));
    return static_cast<uint32_t>(bits >> 32);
}

uint32_t lo_word(const double x) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof(bits));
    return static_cast<uint32_t>(bits);
}

double with_hi(const double x, const uint32_t hi) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof(bits));
    bits = (static_cast<std::uint64_t>(hi) << 32) | (bits & 0xffffffffu);
    double out = 0.0;
    std::memcpy(&out, &bits, sizeof(out));
    return out;
}

// Java's d2i conversion: NaN gives 0, out-of-range values saturate.
int32_t d2i(const double v) {
    if (std::isnan(v)) {
        return 0;
    }
    if (v >= 2147483647.0) {
        return 2147483647;
    }
    if (v <= -2147483648.0) {
        return static_cast<int32_t>(-2147483647 - 1);
    }
    return static_cast<int32_t>(v); // truncates toward zero, as d2i does
}

constexpr double kInf = std::numeric_limits<double>::infinity();

} // namespace

double fdlibm_exp(double x) {
    // FdLibm$Exp static fields.
    static constexpr double half[2] = {0.5, -0.5};
    static constexpr double ln2hi[2] = {0.6931471803691238, -0.6931471803691238};
    static constexpr double ln2lo[2] = {1.9082149292705877e-10, -1.9082149292705877e-10};
    static constexpr double invln2 = 1.4426950408889634;
    static constexpr double p1 = 0.16666666666666602;
    static constexpr double p2 = -0.0027777777777015593;
    static constexpr double p3 = 6.613756321437934e-5;
    static constexpr double p4 = -1.6533902205465252e-6;
    static constexpr double p5 = 4.1381367970572385e-8;
    static constexpr double huge = 1.0e300;
    static constexpr double twom1000 = 9.332636185032189e-302;
    static constexpr double o_threshold = 709.782712893384;
    static constexpr double u_threshold = -745.1332191019411;

    double hi = 0.0;
    double lo = 0.0;
    int32_t n = 0;
    uint32_t hx = hi_word(x);
    const int32_t sgn = static_cast<int32_t>((static_cast<int32_t>(hx) >> 31) & 1);
    uint32_t ix = hx & 0x7fffffffu;
    if (ix >= 0x40862e42u) { // x >= -27*ln2 or NaN/Inf
        if (ix >= 0x7ff00000u) {
            if (((ix & 0xfffffu) | lo_word(x)) != 0) {
                return x + x; // NaN
            }
            return sgn != 0 ? 0.0 : x;
        }
        if (!(x <= o_threshold)) { // dcmpl ifle
            return kInf;
        }
        if (x < u_threshold) { // dcmpg ifge
            return 0.0;
        }
    }
    if (ix > 0x3fd62e42u) { // |x| >= 0.5*ln2
        if (ix < 0x3fe62e42u) { // |x| < ln2
            hi = x - ln2hi[sgn];
            lo = ln2lo[sgn];
            n = 1 - sgn - sgn;
        } else {
            n = d2i(invln2 * x + half[sgn]);
            const double p = static_cast<double>(n);
            hi = x - p * ln2hi[0];
            lo = p * ln2lo[0];
        }
        x = hi - lo;
    } else if (ix >= 0x3e300000u) { // bytecode constant 1043333120
        n = 0;
    } else if (huge + x > 1.0) { // dcmpl ifle
        return 1.0 + x;
    }
    const double z = x * x;
    const double t = x - z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
    if (n == 0) {
        return 1.0 - (x * t / (t - 2.0) - x);
    }
    double y = 1.0 - ((lo - (x * t / (2.0 - t))) - hi);
    if (n >= -1021) {
        // n * 2^20 (|n| <= 1023 stays exact in int): bit-identical to the
        // two's-complement shift the JVM emits, without the C++ UB.
        y = with_hi(y, static_cast<uint32_t>(hi_word(y) + static_cast<uint32_t>(n * 0x100000)));
        return y;
    }
    y = with_hi(y,
        static_cast<uint32_t>(hi_word(y) + static_cast<uint32_t>((n + 1000) * 0x100000)));
    return y * twom1000;
}

double fdlibm_log(double x) {
    // FdLibm$Log static fields.
    static constexpr double ln2_hi = 0.6931471803691238;
    static constexpr double ln2_lo = 1.9082149292705877e-10;
    static constexpr double lg1 = 0.6666666666666735;
    static constexpr double lg2 = 0.3999999999940942;
    static constexpr double lg3 = 0.2857142874366239;
    static constexpr double lg4 = 0.22222198432149784;
    static constexpr double lg5 = 0.1818357216161805;
    static constexpr double lg6 = 0.15313837699209373;
    static constexpr double lg7 = 0.14798198605116586;

    uint32_t hx = hi_word(x);
    const uint32_t lx = lo_word(x);
    int32_t k = 0;
    if (static_cast<int32_t>(hx) < 0x00100000) { // Java compares hx as a signed int
        if ((((hx & 0x7fffffffu)) | lx) != 0) {
            if (static_cast<int32_t>(hx) < 0) {
                return (x - x) / 0.0; // negative argument
            }
            k -= 54;
            x *= 1.8014398509481984e16; // 2^54
            hx = hi_word(x);
        } else {
            return -kInf; // log(0)
        }
    }
    if (hx >= 0x7ff00000u) {
        return x + x; // +Inf or NaN
    }
    k += static_cast<int32_t>(hx >> 20) - 1023;
    int32_t i = static_cast<int32_t>(hx & 0xfffffu);
    int32_t j = static_cast<int32_t>((static_cast<uint32_t>(i) + 614244u) & 0x100000u);
    x = with_hi(x, static_cast<uint32_t>(i) | (static_cast<uint32_t>(j) ^ 0x3ff00000u));
    k += j >> 20;
    const double f = x - 1.0;
    if ((((static_cast<uint32_t>(i) + 2u) & 0xfffffu)) >= 3u) {
        const double s = f / (f + 2.0);
        const double dk = static_cast<double>(k);
        const double z = s * s;
        int32_t i2 = i - 398458;
        const double w = z * z;
        const int32_t i3 = 440401 - i;
        const double r1 = w * (lg2 + w * (lg4 + w * lg6));
        const double r2 = z * (lg1 + w * (lg3 + w * (lg5 + w * lg7)));
        i2 |= i3;
        const double r10 = r2 + r1;
        if (i2 > 0) {
            const double b = (0.5 * f) * f;
            if (k != 0) {
                return dk * ln2_hi - ((b - (s * (b + r10) + dk * ln2_lo)) - f);
            }
            return f - (b - s * (b + r10));
        }
        if (k != 0) {
            return dk * ln2_hi - ((s * (f - r10) - dk * ln2_lo) - f);
        }
        return f - s * (f - r10);
    }
    if (f != 0.0) {
        const double hfsq = (f * f) * (0.5 - (0.3333333333333333 * f));
        if (k != 0) {
            const double dk = static_cast<double>(k);
            return dk * ln2_hi - ((hfsq - dk * ln2_lo) - f);
        }
        return f - hfsq;
    }
    if (k == 0) {
        return 0.0;
    }
    const double dk = static_cast<double>(k);
    return dk * ln2_hi + dk * ln2_lo;
}

// ---------------------------------------------------------------------------
// CI-FDLIBM routing kernels. Provenance: include/kalos/detail/fdlibm.hpp.
// Same operation-for-operation discipline as fdlibm_exp/fdlibm_log above:
// every double expression keeps the operand order of the dumped
// dadd/dsub/dmul/ddiv sequence (the javap dumps under
// tmp_opencode/ci-fdlibm/dumps/ are the source of record; the OpenJDK 21u
// FdLibm.java text was used as a readable cross-check and agrees with them
// everywhere), Java int arithmetic is mirrored on uint32_t with the signed
// comparisons the bytecode makes, and Java d2i semantics reuse the d2i
// helper above.
// ---------------------------------------------------------------------------

double fdlibm_expm1(double x) {
    // java.lang.FdLibm$Expm1 static fields (dump FdExpm1.txt).
    static constexpr double huge = 1.0e300;
    static constexpr double tiny = 1.0e-300;
    static constexpr double o_threshold = 709.7827128933840;
    static constexpr double ln2_hi = 0.6931471803691238;
    static constexpr double ln2_lo = 1.9082149292705877e-10;
    static constexpr double invln2 = 1.4426950408889634;
    static constexpr double q1 = -0.03333333333333313;
    static constexpr double q2 = 0.0015873015872548146;
    static constexpr double q3 = -7.93650757867488e-5;
    static constexpr double q4 = 4.008217827329362e-6;
    static constexpr double q5 = -2.0109921818362437e-7;

    double c = 0.0;
    uint32_t ix = hi_word(x) & 0x7fffffffu;
    const bool negative = (hi_word(x) & 0x80000000u) != 0;
    if (ix >= 0x4043687au) { // |x| >= 56*ln2
        if (ix >= 0x40862e42u) { // |x| >= 709.78...
            if (ix >= 0x7ff00000u) {
                if (((ix & 0xfffffu) | lo_word(x)) != 0) {
                    return x + x; // NaN
                }
                return negative ? -1.0 : x; // exp(+-inf) = {inf, -1}
            }
            if (x > o_threshold) {
                return kInf; // overflow (huge*huge)
            }
        }
        if (negative && (x + tiny) < 0.0) {
            return -1.0; // x < -56*ln2 (tiny - 1.0 in the source)
        }
    }
    int32_t k = 0;
    double hi = 0.0;
    double lo = 0.0;
    if (ix > 0x3fd62e42u) { // |x| > 0.5*ln2
        if (ix < 0x3ff0a2b2u) { // |x| < 1.5*ln2
            if (!negative) {
                hi = x - ln2_hi;
                lo = ln2_lo;
                k = 1;
            } else {
                hi = x + ln2_hi;
                lo = -ln2_lo;
                k = -1;
            }
        } else {
            k = d2i(invln2 * x + (negative ? -0.5 : 0.5)); // signed x, as dumped
            const double t = static_cast<double>(k);
            hi = x - t * ln2_hi; // t*ln2_hi is exact here
            lo = t * ln2_lo;
        }
        x = hi - lo;
        c = (hi - x) - lo;
    } else if (ix < 0x3c900000u) { // |x| < 2^-54
        const double t = huge + x;
        return x - (t - (huge + x));
    }
    // x is now in the primary range
    const double hfx = 0.5 * x;
    const double hxs = x * hfx;
    const double r1 = 1.0 + hxs * (q1 + hxs * (q2 + hxs * (q3 + hxs * (q4 + hxs * q5))));
    const double t = 3.0 - r1 * hfx;
    double e = hxs * ((r1 - t) / (6.0 - x * t));
    if (k == 0) {
        return x - (x * e - hxs); // c is 0
    }
    e = (x * (e - c) - c);
    e -= hxs;
    if (k == -1) {
        return 0.5 * (x - e) - 0.5;
    }
    if (k == 1) {
        if (x < -0.25) {
            return -2.0 * (e - (x + 0.5));
        }
        return 1.0 + 2.0 * (x - e);
    }
    if (k <= -2 || k > 56) { // suffice to return exp(x) - 1
        double y = 1.0 - (e - x);
        y = with_hi(y, static_cast<uint32_t>(hi_word(y)) + (static_cast<uint32_t>(k) << 20));
        return y - 1.0;
    }
    double y = 1.0;
    if (k < 20) {
        y = with_hi(y, 0x3ff00000u - (0x200000u >> static_cast<unsigned>(k))); // y = 1 - 2^-k
        y = y - (e - x);
        y = with_hi(y, static_cast<uint32_t>(hi_word(y)) + (static_cast<uint32_t>(k) << 20));
    } else {
        y = with_hi(y, static_cast<uint32_t>(1023 - k) << 20); // y = 2^-k
        y = x - (e + y);
        y = y + 1.0;
        y = with_hi(y, static_cast<uint32_t>(hi_word(y)) + (static_cast<uint32_t>(k) << 20));
    }
    return y;
}

double fdlibm_tanh(double x) {
    // java.lang.FdLibm$Tanh (dump FdTanh.txt). The pinned kernel is the
    // JDK 9+ expm1-based rewrite, NOT the fdlibm-5.x e_tanh.c; Math.tanh ->
    // StrictMath.tanh -> this kernel, StrictMath.expm1 -> fdlibm_expm1 below.
    const int32_t jx = static_cast<int32_t>(hi_word(x));
    const uint32_t ix = hi_word(x) & 0x7fffffffu;
    if (ix >= 0x7ff00000u) { // x is Inf or NaN
        if (jx >= 0) {
            return 1.0 / x + 1.0;
        }
        return 1.0 / x - 1.0;
    }
    double z = 1.0; // the dumped code stores 1.0 - tiny == 1.0 for |x| >= 22
    if (ix < 0x40360000u) { // |x| < 22
        if (ix < 0x3c800000u) { // |x| < 2^-55
            return x * (1.0 + x);
        }
        if (ix >= 0x3ff00000u) { // |x| >= 1
            const double t = fdlibm_expm1(2.0 * std::fabs(x));
            z = 1.0 - 2.0 / (t + 2.0);
        } else {
            const double t = fdlibm_expm1(-2.0 * std::fabs(x));
            z = -t / (t + 2.0);
        }
    }
    return jx >= 0 ? z : -z;
}

double fdlibm_log1p(double x) {
    // java.lang.FdLibm$Log1p static fields (dump FdLog1p.txt).
    static constexpr double ln2_hi = 0.6931471803691238;
    static constexpr double ln2_lo = 1.9082149292705877e-10;
    static constexpr double lp1 = 0.6666666666666735;
    static constexpr double lp2 = 0.3999999999940942;
    static constexpr double lp3 = 0.2857142874366239;
    static constexpr double lp4 = 0.22222198432149784;
    static constexpr double lp5 = 0.1818357216161805;
    static constexpr double lp6 = 0.15313837699209373;
    static constexpr double lp7 = 0.14798198605116586;
    static constexpr double two54 = 1.8014398509481984e16;

    double f = 0.0;
    double c = 0.0;
    int32_t hu = 0;
    const int32_t hx = static_cast<int32_t>(hi_word(x));
    const uint32_t ax = hi_word(x) & 0x7fffffffu;
    int32_t k = 1;
    if (hx < 1071284858) { // hx < 0x3FDA827A: x < 0.41422 (all negatives included)
        if (ax >= 0x3ff00000u) { // x <= -1.0
            if (x == -1.0) {
                return -kInf; // log1p(-1) = -inf
            }
            return std::numeric_limits<double>::quiet_NaN(); // log1p(x < -1)
        }
        if (ax < 0x3e200000u) { // |x| < 2^-29
            if (two54 + x > 0.0 && ax < 0x3c900000u) { // |x| < 2^-54
                return x;
            }
            return x - x * x * 0.5;
        }
        if (hx > 0 || hx <= -1076707645) { // 0xBFD2BEC3: -0.2929 < x < 0.41422
            k = 0;
            f = x;
            hu = 1;
        }
    }
    if (hx >= 0x7ff00000) {
        return x + x; // +Inf or NaN
    }
    if (k != 0) {
        double u = 0.0;
        if (hx < 1128267776) { // hx < 0x43400000
            u = 1.0 + x;
            hu = static_cast<int32_t>(hi_word(u));
            k = (hu >> 20) - 1023;
            c = (k > 0) ? 1.0 - (u - x) : x - (u - 1.0); // correction term
            c /= u;
        } else {
            u = x;
            hu = static_cast<int32_t>(hi_word(u));
            k = (hu >> 20) - 1023;
            c = 0.0;
        }
        hu &= 0x000fffff;
        if (hu < 434334) { // 0x6A09E
            u = with_hi(u, static_cast<uint32_t>(hu) | 0x3ff00000u); // normalize u
        } else {
            k += 1;
            u = with_hi(u, static_cast<uint32_t>(hu) | 0x3fe00000u); // normalize u/2
            hu = (1048576 - hu) >> 2;
        }
        f = u - 1.0;
    }
    const double hfsq = 0.5 * f * f; // (0.5*f)*f, as dumped
    if (hu == 0) { // |f| < 2^-20
        if (f == 0.0) {
            if (k == 0) {
                return 0.0;
            }
            c = c + static_cast<double>(k) * ln2_lo;
            return static_cast<double>(k) * ln2_hi + c;
        }
        const double r = hfsq * (1.0 - 0.6666666666666666 * f);
        if (k == 0) {
            return f - r;
        }
        return static_cast<double>(k) * ln2_hi - ((r - (static_cast<double>(k) * ln2_lo + c)) - f);
    }
    const double s = f / (2.0 + f);
    const double z = s * s;
    const double r = z * (lp1 + z * (lp2 + z * (lp3 + z * (lp4 + z * (lp5 + z * (lp6 + z * lp7))))));
    if (k == 0) {
        return f - (hfsq - s * (hfsq + r));
    }
    return static_cast<double>(k) * ln2_hi -
        ((hfsq - (s * (hfsq + r) + (static_cast<double>(k) * ln2_lo + c))) - f);
}

double fdlibm_atanh(double x) {
    // fdlibm-5.x __ieee754_atanh (netlib.org/fdlibm/e_atanh.c 1.3 95/01/18,
    // fetched 2026-11-21), the lineage sibling of the vendored kernels. NO
    // JVM counterpart exists: neither java.lang.Math nor StrictMath has atanh
    // in the pinned JDK 21 (verified by live reflection + javap). This kernel
    // pins the reflection-coordinate map (arctanh -> theta) to
    // machine-independent bits; it carries no bit-exactness claim against the
    // JVM. See the header provenance table.
    static constexpr double one = 1.0;
    static constexpr double huge = 1.0e300;
    static constexpr double zero = 0.0;
    const int32_t hx = static_cast<int32_t>(hi_word(x));
    const uint32_t lx = lo_word(x);
    const uint32_t ix = hi_word(x) & 0x7fffffffu;
    if ((ix | ((lx | (0u - lx)) >> 31)) > 0x3ff00000u) { // |x| > 1
        return (x - x) / (x - x); // NaN
    }
    if (ix == 0x3ff00000u) { // |x| == 1 (exact: lx zero, see the test above)
        return x / zero; // atanh(+-1) = +-inf
    }
    if (ix < 0x3e300000u && (huge + x) > zero) { // |x| < 2^-27
        return x;
    }
    const double y = with_hi(x, ix); // |x|
    double t = 0.0;
    if (ix < 0x3fe00000u) { // |x| < 1/2
        t = y + y;
        t = 0.5 * fdlibm_log1p(t + t * y / (one - y));
    } else {
        t = 0.5 * fdlibm_log1p((y + y) / (one - y));
    }
    return hx >= 0 ? t : -t;
}

namespace {

// java.lang.FdLibm$RemPio2 / $Sin / $Cos constant set (dumps
// FdRemPio2.txt, FdSin.txt, FdCos.txt).
constexpr double kInvpio2 = 0.6366197723675814;
constexpr double kPio2_1 = 1.5707963267341256;
constexpr double kPio2_1t = 6.077100506506192e-11;
constexpr double kPio2_2 = 6.077100506303966e-11;
constexpr double kPio2_2t = 2.0222662487959506e-21;
constexpr double kPio2_3 = 2.0222662487111665e-21;
constexpr double kPio2_3t = 8.4784276603689e-32;

// FdLibm$RemPio2 npio2_hw (int[32], extracted from the static initializer of
// the dump; verified bit-identical to the OpenJDK 21u table).
const uint32_t kNpio2_hw[32] = {
    0x3ff921fbu, 0x400921fbu, 0x4012d97cu, 0x401921fbu, 0x401f6a7au, 0x4022d97cu,
    0x4025fdbbu, 0x402921fbu, 0x402c463au, 0x402f6a7au, 0x4031475cu, 0x4032d97cu,
    0x40346b9cu, 0x4035fdbbu, 0x40378fdbu, 0x403921fbu, 0x403ab41bu, 0x403c463au,
    0x403dd85au, 0x403f6a7au, 0x40407e4cu, 0x4041475cu, 0x4042106cu, 0x4042d97cu,
    0x4043a28cu, 0x40446b9cu, 0x404534acu, 0x4045fdbbu, 0x4046c6cbu, 0x40478fdbu,
    0x404858ebu, 0x404921fbu,
};

double kernel_cos(double x, double y) {
    // FdLibm$Cos.__kernel_cos (dump FdCos.txt).
    static constexpr double c1 = 0.0416666666666666;
    static constexpr double c2 = -0.001388888888887411;
    static constexpr double c3 = 2.480158728947673e-5;
    static constexpr double c4 = -2.7557314351390663e-7;
    static constexpr double c5 = 2.087572321298175e-9;
    static constexpr double c6 = -1.1359647557788195e-11;
    const uint32_t ix = hi_word(x) & 0x7fffffffu;
    if (ix < 0x3e400000u && d2i(x) == 0) {
        return 1.0;
    }
    const double z = x * x;
    const double r = z * (c1 + z * (c2 + z * (c3 + z * (c4 + z * (c5 + z * c6)))));
    if (ix < 0x3fd33333u) { // |x| < 0.3
        return 1.0 - (0.5 * z - (z * r - x * y));
    }
    double qx = 0.0;
    if (ix > 0x3fe90000u) { // x > 0.78125
        qx = 0.28125;
    } else {
        qx = with_hi(0.0, ix - 0x00200000u); // __HI_LO(ix - 0x200000, 0)
    }
    const double hz = 0.5 * z - qx;
    const double a = 1.0 - qx;
    return a - (hz - (z * r - x * y));
}

double kernel_sin(double x, double y, int iy) {
    // FdLibm$Sin.__kernel_sin (dump FdSin.txt).
    static constexpr double s1 = -0.16666666666666632;
    static constexpr double s2 = 0.00833333333332249;
    static constexpr double s3 = -1.984126982985795e-4;
    static constexpr double s4 = 2.7557313707070068e-6;
    static constexpr double s5 = -2.5050760253406863e-8;
    static constexpr double s6 = 1.58969099521155e-10;
    const uint32_t ix = hi_word(x) & 0x7fffffffu;
    if (ix < 0x3e400000u && d2i(x) == 0) {
        return x;
    }
    const double z = x * x;
    const double v = z * x;
    const double r = s2 + z * (s3 + z * (s4 + z * (s5 + z * s6)));
    if (iy == 0) {
        return x + v * (s1 + z * r);
    }
    return x - ((z * (0.5 * y - v * r) - y) - v * s1);
}

// FdLibm$RemPio2.__ieee754_rem_pio2 (dump FdRemPio2.txt), complete through
// the medium-size Cody-Waite branch (|x| high word <= 0x413921FB). The
// huge-argument branch dispatches to FdLibm$KernelRemPio2 (Payne-Hanek),
// which is NOT transcribed: every routed call site is an angle bounded far
// below 0x1.921fbp+20 (~1.65e6) — UnitRootsSolver.pnorm passes theta in
// (0, 2*pi] and the rational-filter evaluation passes |x| <= pi — so the
// fail-closed throw below is unreachable from the qualified numerics.
int rem_pio2(double x, double* y) {
    const int32_t hx = static_cast<int32_t>(hi_word(x));
    const uint32_t ix = hi_word(x) & 0x7fffffffu;
    if (ix <= 0x3fe921fbu) { // |x| <= pi/4: no reduction
        y[0] = x;
        y[1] = 0.0;
        return 0;
    }
    if (ix < 0x4002d97cu) { // |x| < 3pi/4, special case n = +-1
        if (hx > 0) {
            double z = x - kPio2_1;
            if (ix != kNpio2_hw[0]) { // 33+53 bit pi is good enough
                y[0] = z - kPio2_1t;
                y[1] = (z - y[0]) - kPio2_1t;
            } else { // near pi/2, use 33+33+53 bit pi
                z -= kPio2_2;
                y[0] = z - kPio2_2t;
                y[1] = (z - y[0]) - kPio2_2t;
            }
            return 1;
        }
        double z = x + kPio2_1;
        if (ix != kNpio2_hw[0]) {
            y[0] = z + kPio2_1t;
            y[1] = (z - y[0]) + kPio2_1t;
        } else {
            z += kPio2_2;
            y[0] = z + kPio2_2t;
            y[1] = (z - y[0]) + kPio2_2t;
        }
        return -1;
    }
    if (ix <= 0x413921fbu) { // medium size, Cody-Waite with staged cancellation fixup
        double t = std::fabs(x);
        int n = d2i(t * kInvpio2 + 0.5);
        const double fn = static_cast<double>(n);
        double r = t - fn * kPio2_1;
        double w = fn * kPio2_1t; // 1st round good to 85 bit
        if (n < 32 && ix != kNpio2_hw[n - 1]) { // quick check: no cancellation
            y[0] = r - w;
        } else {
            const uint32_t j = ix >> 20;
            y[0] = r - w;
            int32_t i = static_cast<int32_t>(j) -
                (static_cast<int32_t>(hi_word(y[0])) >> 20 & 0x7ff);
            if (i > 16) { // 2nd iteration needed, good to 118 bit
                t = r;
                w = fn * kPio2_2;
                r = t - w;
                w = fn * kPio2_2t - ((t - r) - w);
                y[0] = r - w;
                i = static_cast<int32_t>(j) -
                    (static_cast<int32_t>(hi_word(y[0])) >> 20 & 0x7ff);
                if (i > 49) { // 3rd iteration needed, 151 bit accuracy
                    t = r; // will cover all possible cases
                    w = fn * kPio2_3;
                    r = t - w;
                    w = fn * kPio2_3t - ((t - r) - w);
                    y[0] = r - w;
                }
            }
        }
        y[1] = (r - y[0]) - w;
        if (hx < 0) {
            y[0] = -y[0];
            y[1] = -y[1];
            return -n;
        }
        return n;
    }
    if (ix >= 0x7ff00000u) { // Inf or NaN: y = NaN (never reached from sin/cos)
        y[0] = x - x;
        y[1] = x - x;
        return 0;
    }
    throw std::domain_error(
        "fdlibm sin/cos: Payne-Hanek argument reduction (|x| > 1.64e6) is not "
        "transcribed; every routed call site is angle-bounded");
}

} // namespace

double fdlibm_cos(double x) {
    // java.lang.FdLibm$Cos.compute (dump FdCos.txt).
    const uint32_t ix = hi_word(x) & 0x7fffffffu;
    if (ix <= 0x3fe921fbu) { // |x| <= pi/4
        return kernel_cos(x, 0.0);
    }
    if (ix >= 0x7ff00000u) { // cos(Inf or NaN) is NaN
        return x - x;
    }
    double y[2] = {0.0, 0.0};
    const int n = rem_pio2(x, y);
    switch (n & 3) {
    case 0:
        return kernel_cos(y[0], y[1]);
    case 1:
        return -kernel_sin(y[0], y[1], 1);
    case 2:
        return -kernel_cos(y[0], y[1]);
    default:
        return kernel_sin(y[0], y[1], 1);
    }
}

double fdlibm_sin(double x) {
    // java.lang.FdLibm$Sin.compute (dump FdSin.txt).
    const uint32_t ix = hi_word(x) & 0x7fffffffu;
    if (ix <= 0x3fe921fbu) { // |x| <= pi/4
        return kernel_sin(x, 0.0, 0);
    }
    if (ix >= 0x7ff00000u) { // sin(Inf or NaN) is NaN
        return x - x;
    }
    double y[2] = {0.0, 0.0};
    const int n = rem_pio2(x, y);
    switch (n & 3) {
    case 0:
        return kernel_sin(y[0], y[1], 1);
    case 1:
        return kernel_cos(y[0], y[1]);
    case 2:
        return -kernel_sin(y[0], y[1], 1);
    default:
        return -kernel_cos(y[0], y[1]);
    }
}


// FdLibm.Pow for the binary64 exponent 1/3 used by DoornikHansen.
// The split logarithm and exponential follow fdlibm e_pow.c:
// Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
// Permission to use, copy, modify, and distribute this software is freely
// granted, provided that this notice is preserved.
double fdlibm_power_one_third(double x) {
    if (x == 0. || x == 1. || std::isinf(x)) return std::abs(x);
    if (!(x > 0.)) return std::numeric_limits<double>::quiet_NaN();
    auto upper = [](double v) { return with_hi(0., hi_word(v)); };

    // Range reduction about 1 or 1.5, with log2(1.5) split at 24 bits.
    int exponent = 0;
    if (hi_word(x) < 0x00100000u) { x *= 0x1p53; exponent -= 53; }
    const auto high = hi_word(x);
    exponent += static_cast<int>(high >> 20) - 1023;
    const auto fraction = high & 0xfffffu;
    auto normalized = fraction | 0x3ff00000u;
    const bool center_one_half = fraction > 0x3988eu && fraction < 0xbb67au;
    if (fraction >= 0xbb67au) { ++exponent; normalized -= 0x100000u; }
    x = with_hi(x, normalized);
    const double center = center_one_half ? 1.5 : 1.;
    const double offset_hi = center_one_half ? 0x1.2b8034p-1 : 0.;
    const double offset_lo = center_one_half ? 0x1.cfdeb43cfd006p-27 : 0.;
    const double diff = x - center, inverse_sum = 1. / (x + center);
    const double ratio = diff * inverse_sum, ratio_hi = upper(ratio);
    const double sum_hi = with_hi(0., ((normalized >> 1) | 0x20000000u) + 0x80000u +
                                     (center_one_half ? 0x40000u : 0u));
    const double sum_lo = x - (sum_hi - center);
    const double ratio_lo = inverse_sum * ((diff - ratio_hi * sum_hi) - ratio_hi * sum_lo);
    const double sq = ratio * ratio;
    double tail = sq * sq * (0x1.3333333333303p-1 + sq * (0x1.b6db6db6fabffp-2 +
        sq * (0x1.55555518f264dp-2 + sq * (0x1.17460a91d4101p-2 +
        sq * (0x1.d864a93c9db65p-3 + sq * 0x1.a7e284a454eefp-3)))));
    tail += ratio_lo * (ratio_hi + ratio);
    const double sq_hi = ratio_hi * ratio_hi;
    const double poly_hi = upper(3. + sq_hi + tail);
    const double poly_lo = tail - ((poly_hi - 3.) - sq_hi);
    const double u = ratio_hi * poly_hi, v = ratio_lo * poly_hi + poly_lo * ratio;
    const double product_hi = upper(u + v), product_lo = v - (product_hi - u);
    const double log_part_hi = 0x1.ec709ep-1 * product_hi;
    const double log_part_lo = -0x1.e2fe0145b01f5p-28 * product_hi +
        product_lo * 0x1.ec709dc3a03fdp-1 + offset_lo;
    const double log_hi = upper(((log_part_hi + log_part_lo) + offset_hi) + exponent);
    const double log_lo = log_part_lo - (((log_hi - exponent) - offset_hi) - log_part_hi);

    constexpr double power = 1. / 3.;
    const double power_hi = upper(power);
    const double exp_lo = (power - power_hi) * log_hi + power * log_lo;
    double exp_hi = power_hi * log_hi;
    const double argument = exp_lo + exp_hi;
    // FdLibm rounds the power of two once, before the final exponential.
    const auto abs_hi = hi_word(argument) & 0x7fffffffu;
    int scale = 0;
    if (abs_hi > 0x3fe00000u) {
        const int k = static_cast<int>(abs_hi >> 20) - 1023;
        const auto rounded = abs_hi + (0x100000u >> (k + 1));
        const int width = static_cast<int>(rounded >> 20) - 1023;
        scale = static_cast<int>(((rounded & 0xfffffu) | 0x100000u) >> (20 - width));
        if (argument < 0.) scale = -scale;
        exp_hi -= scale;
    }
    const double t = upper(exp_lo + exp_hi);
    const double z_hi = t * 0x1.62e43p-1;
    const double z_lo = (exp_lo - (t - exp_hi)) * 0x1.62e42fefa39efp-1 - t * 0x1.05c610ca86c39p-29;
    const double z = z_hi + z_lo, w = z_lo - (z - z_hi), z2 = z * z;
    const double a = z - z2 * (0x1.555555555553ep-3 + z2 * (-0x1.6c16c16bebd93p-9 +
        z2 * (0x1.1566aaf25de2cp-14 + z2 * (-0x1.bbd41c5d26bf1p-20 + z2 * 0x1.6376972bea4d0p-25))));
    const double remainder = (z * a) / (a - 2.) - (w + z * w);
    return std::scalbn(1. - (remainder - z), scale);
}

} // namespace kalos::detail
