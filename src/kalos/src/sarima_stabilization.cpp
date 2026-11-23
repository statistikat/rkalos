#include "kalos/detail/sarima_stabilization.hpp"

#include "kalos/linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace kalos::detail {

// FastEigenValuesSolver uses the closed quadratic formula through degree two.
std::vector<ComplexRoot> jvm_closed_roots(const std::vector<double>& asc) {
    if (asc.empty()) return {};
    int deg = static_cast<int>(asc.size()) - 1;
    while (deg > 0 && asc[deg] == 0.0) {
        --deg;
    }
    std::vector<ComplexRoot> out;
    if (deg == 0) {
        return out;
    }
    if (deg == 1) {
        out.push_back({-asc[0] / asc[1], 0.0});
        return out;
    }
    if (deg > 2) {
        throw std::invalid_argument("closed-form roots require degree at most two");
    }
    const double a = asc[2], b = asc[1], c0 = asc[0];
    const double denom = 2.0 * a;
    const double disc = b * b - 4.0 * a * c0;
    if (disc < 0.0) {
        const double s = std::sqrt(-disc);
        const ComplexRoot z{-b / denom, s / denom};
        out.push_back(z);
        out.push_back({z.real, -z.imag});
        return out;
    }
    const double s = std::sqrt(disc);
    out.push_back({(-b + s) / denom, 0.0});
    out.push_back({(-b - s) / denom, 0.0});
    return out;
}

double jvm_complex_abs(double re, double im) {
    if (re == 0.0 && im == 0.0) {
        return 0.0;
    }
    const double a = std::fabs(re), b = std::fabs(im);
    const double m = std::max(a, b), n = std::min(a, b);
    if (n == 0.0) {
        return m;
    }
    const double t = n / m;
    return m * std::sqrt(1.0 + t * t);
}

bool stabilize_coefficient_block(std::vector<double>& block) {
    std::size_t length = block.size();
    while (length > 0 && block[length - 1] == 0.0) --length;
    if (length == 0) return false;
    if (length == 1) {
        if (std::abs(block[0]) > 1.0) {
            block[0] = 1.0 / block[0];
            return true;
        }
        return false;
    }
    std::vector<double> polynomial(length + 1, 0.0);
    polynomial[0] = 1.0;
    for (std::size_t i = 0; i < length; ++i) polynomial[i + 1] = block[i];
    const auto found = length <= 2 ? jvm_closed_roots(polynomial) : polynomial_roots(polynomial);
    if (found.empty()) return false;
    std::vector<std::complex<double>> roots;
    roots.reserve(found.size());
    for (const auto& root : found) roots.emplace_back(root.real, root.imag);
    bool changed = false;
    for (auto& root : roots) {
        const double re = root.real(), im = root.imag();
        const double modulus = jvm_complex_abs(re, im);
        if (1.0 / modulus > 1.0) {
            if (modulus == 0.0) {
                throw std::runtime_error("stabilization: zero root cannot be inverted");
            }
            if (std::abs(re) >= std::abs(im)) {
                const double inv = 1.0 / (re + im * (im / re));
                root = {inv, inv * (-im / re)};
            } else {
                const double inv = 1.0 / (re * (re / im) + im);
                root = {inv * (re / im), -inv};
            }
            changed = true;
        }
    }
    if (!changed) return false;
    // Pinned Polynomial.fromComplexRoots followed by division by the
    // constant term; only real parts cross back into the block.
    std::vector<std::complex<double>> z(roots.size() + 1, {0.0, 0.0});
    z[0] = {1.0, 0.0};
    for (std::size_t i = 0; i < roots.size(); ++i) {
        z[i + 1] = z[i];
        for (std::size_t j = i; j >= 1; --j) z[j] = z[j - 1] - z[j] * roots[i];
        z[0] = z[0] * (-roots[i]);
    }
    if (!(z[0].real() != 0.0)) {
        throw std::runtime_error("stabilization: roots must all be different from 0");
    }
    for (std::size_t i = 0; i < length; ++i) {
        if (!std::isfinite(z[i + 1].real())) {
            throw std::runtime_error("stabilization: unstable reconstruction");
        }
        block[i] = z[i + 1].real() * (1.0 / z[0].real());
    }
    return true;
}

} // namespace kalos::detail
