// Item-29 probe slice: deterministic synthetic-series generator for the
// frozen corpora. Uses a fixed 64-bit LCG (Knuth MMIX + mix64 output) plus
// polar-free Box-Muller so emitted CSV bytes are reproducible on any
// IEEE-754 binary64 / two's-complement uint64 platform — no host RNG
// version sensitivity. Values printed with %.17g.
//
// Model (kalos convention):
//   a(B) x_t = c(B) e_t,
//   a(B) = (1 - sum ar_i B^i)(1 - sum sar_j B^(s j)),
//   c(B) = (1 + sum ma_i B^i)(1 + sum sma_j B^(s j)),  e_t ~ N(0,1) iid
// then y = (1-B)^(-d) (1-B^s)^(-D) x; the ARMA warm-up runs 4000 burn-in
// steps that are discarded.
//
// usage: gen_series out.csv freq d D "ar,ar2;ma,ma2;sar;sma" n seed
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace {

struct Lcg {
    std::uint64_t s;
    explicit Lcg(std::uint64_t seed) : s(seed) {}
    std::uint64_t next() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        std::uint64_t z = s;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    double uniform() { // dyadic in (0,1)
        return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

std::vector<double> parse_vec(const std::string& s) {
    std::vector<double> v;
    std::string item;
    for (std::size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == ',') {
            if (!item.empty()) v.push_back(std::stod(item));
            item.clear();
        } else item.push_back(s[i]);
    }
    return v;
}

std::vector<double> mul_factor(const std::vector<double>& p, const std::vector<double>& coef,
                               int step, int sign) {
    // multiply polynomial p(B) by prod_j (1 + sign * coef_j-1 B^(step j))
    std::vector<double> q(p.size() + static_cast<std::size_t>(step) * coef.size(), 0.0);
    for (std::size_t i = 0; i < p.size(); ++i) q[i] += p[i];
    for (std::size_t j = 0; j < coef.size(); ++j)
        for (std::size_t i = 0; i < p.size(); ++i)
            q[i + static_cast<std::size_t>(step) * (j + 1)] += sign * coef[j] * p[i];
    return q;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 8) {
        std::fprintf(stderr, "usage: gen_series out.csv freq d D \"ar;ma;sar;sma\" n seed\n");
        return 2;
    }
    const std::string out_path = argv[1];
    const int freq = std::atoi(argv[2]);
    const int d = std::atoi(argv[3]);
    const int D = std::atoi(argv[4]);
    std::vector<std::string> tokens;
    {
        std::string cur;
        for (const char c : std::string(argv[5])) {
            if (c == ';') { tokens.push_back(cur); cur.clear(); }
            else cur.push_back(c);
        }
        tokens.push_back(cur);
    }
    std::vector<double> ar_v, ma_v, sar_v, sma_v;
    if (tokens.size() > 0) ar_v = parse_vec(tokens[0]);
    if (tokens.size() > 1) ma_v = parse_vec(tokens[1]);
    if (tokens.size() > 2) sar_v = parse_vec(tokens[2]);
    if (tokens.size() > 3) sma_v = parse_vec(tokens[3]);
    const int n = std::atoi(argv[6]);
    const std::uint64_t seed = std::strtoull(argv[7], nullptr, 10);

    std::vector<double> a = {1.0}, c = {1.0};
    a = mul_factor(a, ar_v, 1, -1);
    a = mul_factor(a, sar_v, freq, -1);
    c = mul_factor(c, ma_v, 1, +1);
    c = mul_factor(c, sma_v, freq, +1);
    const int la = static_cast<int>(a.size()) - 1;
    const int lc = static_cast<int>(c.size()) - 1;
    const int burn = 4000;
    const int total = burn + n;

    Lcg rng(seed);
    std::vector<double> x(static_cast<std::size_t>(total), 0.0);
    std::vector<double> e(static_cast<std::size_t>(total), 0.0);
    for (int t = 0; t < total; ++t) {
        const double u1 = rng.uniform();
        const double u2 = rng.uniform();
        const double noise = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        e[static_cast<std::size_t>(t)] = noise;
        double value = noise;
        for (int j = 1; j <= la; ++j)
            if (t >= j) value += a[j] * x[static_cast<std::size_t>(t - j)];
        for (int j = 1; j <= lc; ++j)
            if (t >= j) value += c[j] * e[static_cast<std::size_t>(t - j)];
        x[static_cast<std::size_t>(t)] = value;
    }

    std::vector<double> y;
    for (int t = burn; t < total; ++t) y.push_back(x[static_cast<std::size_t>(t)]);
    auto integrate = [&](int lag) {
        std::vector<double> z = y;
        for (int t = lag; t < n; ++t) z[t] += z[t - lag];
        y = z;
    };
    for (int i = 0; i < d; ++i) integrate(1);
    for (int i = 0; i < D; ++i) integrate(freq);

    std::ofstream out(out_path);
    out << "y\n";
    for (const auto value : y) {
        char buffer[64];
        std::snprintf(buffer, sizeof buffer, "%.17g", value);
        out << buffer << "\n";
    }
    out.close();
    return 0;
}
