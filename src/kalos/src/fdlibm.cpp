#include "kalos/detail/fdlibm.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

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
        y = with_hi(y, static_cast<uint32_t>(hi_word(y) + static_cast<uint32_t>(n << 20)));
        return y;
    }
    y = with_hi(y,
        static_cast<uint32_t>(hi_word(y) + static_cast<uint32_t>((n + 1000) << 20)));
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

} // namespace kalos::detail
