#include "kalos/calendar_generation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace kalos {
namespace {

// ---------------------------------------------------------------- domains --

int ratio_of(int frequency) { return 12 / frequency; }  // months per period

void validate_domain(const calendar_domain& d) {
    if (d.frequency != 12 && d.frequency != 4 && d.frequency != 1)
        throw std::invalid_argument(
            "calendar_generation: domain frequency must be 12, 4 or 1");
    if (d.length <= 0)
        throw std::invalid_argument(
            "calendar_generation: domain length must be positive");
    if (d.year < 1 || d.year > 9999)
        throw std::invalid_argument(
            "calendar_generation: domain year out of range [1, 9999]");
    int m = d.month;
    if (d.frequency == 12) {
        if (m < 1 || m > 12)
            throw std::invalid_argument(
                "calendar_generation: domain month not aligned with frequency");
    } else if (d.frequency == 4) {
        if (m != 1 && m != 4 && m != 7 && m != 10)
            throw std::invalid_argument(
                "calendar_generation: domain month not aligned with frequency");
    } else if (m != 1) {
        throw std::invalid_argument(
            "calendar_generation: domain month not aligned with frequency");
    }
}

void require_trading_frequency(const calendar_domain& d, const char* what) {
    if (d.frequency != 12 && d.frequency != 4)
        throw std::invalid_argument(std::string("calendar_generation: ") + what +
                                    " requires a monthly or quarterly domain");
}

// floor division (Java TsPeriod positions floor toward -inf)
int fdiv(int a, int b) {
    int q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
    return q;
}

// period index of calendar month (y, m) relative to the domain, or a
// refusal-free position (may be < 0 or >= length; callers guard).
int pos_of(const calendar_domain& d, int y, int m) {
    int dm = (y * 12 + (m - 1)) - (d.year * 12 + (d.month - 1));
    return fdiv(dm, ratio_of(d.frequency));
}

bool is_leap(int y) { return y % 4 == 0 && (y % 100 != 0 || y % 400 == 0); }

int month_len(int y, int m) {
    static const int L[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return L[m - 1] + (m == 2 && is_leap(y) ? 1 : 0);
}

// Hinnant days_from_civil (epoch 1970-01-01) and inverse.
long days_from_civil(int y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return static_cast<long>(era) * 146097 + static_cast<long>(doe) - 719468;
}

// Monday-based weekday (0=Mon..6=Sun); 1970-01-01 was a Thursday (3).
int weekday_mon(int y, int m, int d) {
    long z = days_from_civil(y, static_cast<unsigned>(m), static_cast<unsigned>(d));
    long w = (z + 3) % 7;  // +3: epoch day 0 = Thursday = index 3
    if (w < 0) w += 7;
    return static_cast<int>(w);
}

// --------------------------------------------------------------- tables ---

const double MDAYS[12] = {31.0, 28.25, 31.0, 30.0, 31.0, 30.0,
                          31.0, 31.0,   30.0, 31.0, 30.0, 31.0};

const double MEANS_MAR[25] = {
    0.6576,     0.6450833, 0.6311304, 0.6162727, 0.5999048, 0.583,
    0.5661053,  0.549,     0.5318824, 0.514625,  0.4973333, 0.4807143,
    0.4643077,  0.4476667, 0.4305455, 0.4136,    0.3975556, 0.382,
    0.3654286,  0.3483333, 0.3304,    0.3125,    0.2966667, 0.281,
    0.266};
const double MEANS_APR[25] = {
    0.33872,    0.3528333, 0.3677391, 0.3834545, 0.4000952, 0.417,
    0.4338947,  0.451,     0.4681176, 0.485375,  0.5026667, 0.5192857,
    0.5356923,  0.5523333, 0.5694545, 0.5864,    0.6024444, 0.618,
    0.6345714,  0.6516667, 0.6696,    0.6875,    0.7033333, 0.719,
    0.734};

const int C_MAR2[28] = {1,    3,    8,    22,   51,   91,   145,  226,  343,
                        488,  657,  858,  1099, 1379, 1692, 2036, 2414, 2828,
                        3284, 3779, 4307, 4866, 5458, 6092, 6764, 7469, 8204,
                        8971};
const int C_APR2[28] = {878,   1784,  2713,  3664,  4633,  5613,  6586,  7547,
                        8484,  9400,  10293, 11155, 11978, 12762, 13513, 14233,
                        14919, 15569, 16177, 16746, 17282, 17787, 18259, 18689,
                        19081, 19440, 19769, 20066};
const int C_MAY2[28] = {185, 341, 471, 570, 636, 680, 717, 739, 749, 752, 754, 755,
                        755, 755, 755, 755, 755, 755, 755, 755, 755, 755, 755, 755,
                        755, 755, 755, 755};

const int JD[43] = {1,  1,  3,  9,  15, 11, 14, 27, 36, 28, 24, 32, 40, 39,
                    33, 31, 34, 36, 42, 39, 33, 31, 33, 42, 38, 33, 30, 32,
                    39, 40, 33, 29, 26, 31, 33, 22, 7,  15, 12, 7,  1,  1,
                    1};

// DayClustering groupsDefinition[7] (Mon..Sun -> group id), ascending-id
// order defines the contrast column order; the reference is group 0.
struct cluster_info {
    day_clustering clustering;
    const char* tag;
    const char* name;
    int groups;
    int def[7];
};
const cluster_info CLUSTERS[9] = {
    {day_clustering::td2, "TD2", "td2", 2, {1, 1, 1, 1, 1, 0, 0}},
    {day_clustering::td2c, "TD2c", "td2c", 2, {1, 1, 1, 1, 1, 1, 0}},
    {day_clustering::td2d, "TDuser", "tduser", 2, {1, 1, 1, 1, 0, 0, 0}},
    {day_clustering::td3, "TD3", "td3", 3, {1, 1, 1, 1, 1, 2, 0}},
    {day_clustering::td3c, "TD3c", "td3c", 3, {1, 1, 1, 1, 2, 2, 0}},
    {day_clustering::td4, "TD4", "td4", 4, {1, 1, 1, 1, 2, 3, 0}},
    {day_clustering::td4c, "TD4c", "td4c", 4, {1, 2, 2, 2, 2, 3, 0}},
    {day_clustering::td6, "TD6", "td6", 6, {1, 2, 3, 4, 5, 0, 0}},
    {day_clustering::td7, "TD7", "td7", 7, {1, 2, 3, 4, 5, 6, 0}},
};

const cluster_info& cluster_at(day_clustering clustering) {
    for (const auto& c : CLUSTERS)
        if (c.clustering == clustering) return c;
    throw std::invalid_argument("calendar_generation: unknown day clustering");
}

// -------------------------------------------------------------- computus --

}  // namespace

// Easter.easter (Meeus/Jones/Butcher; the factory's only Gregorian computus,
// easter == easter2 pinned 1985..2015).
calendar_date easter_date(int year) {
    if (year < 1 || year > 9999)
        throw std::invalid_argument(
            "calendar_generation: computus year out of range [1, 9999]");
    int a = year % 19;
    int b = year / 100;
    int c = year % 100;
    int g = b / 4;
    int i = b % 4;
    int h = (b + 8) / 25;
    int k = (b - h + 1) / 3;
    int l = (19 * a + b - g - k + 15) % 30;
    int m = c / 4;
    int n = c % 4;
    int p = (32 + 2 * i + 2 * m - l - n) % 7;
    int q = (a + 11 * l + 22 * p) / 451;
    int s = l + p - 7 * q + 114;
    return {year, s / 31, s % 31 + 1};
}

namespace {
// Easter.julianDate2JDN (Julian calendar date -> JDN).
int julian_date_to_jdn(int year, int month, int day) {
    int k = (14 - month) / 12;
    int yy = year + 4800 - k;
    int mm = month + 12 * k - 3;
    return day + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - 32083;
}

// Easter.JDN2GregorianDate (Fliegel-Van Flandern; Java integer division).
calendar_date jdn_to_gregorian(int jd) {
    int j = jd + 1401 + (((4 * jd + 274277) / 146097) * 3) / 4 - 38;
    int m = 4 * j + 3;
    int n = (m % 1461) / 4;
    int i = 5 * n + 2;
    int d = (i % 153) / 5;
    int mon = (i / 153 + 2) % 12;
    int y = m / 1461 - 4716 + (14 - mon) / 12;
    return {y, mon + 1, d + 1};
}
}  // namespace

calendar_date julian_easter_date(int year, bool gregorian_dates) {
    if (year < 1 || year > 9999)
        throw std::invalid_argument(
            "calendar_generation: computus year out of range [1, 9999]");
    int a = year % 19;
    int b = year % 4;
    int c = year % 7;
    int d = (19 * a + 15) % 30;
    int e = (2 * b + 4 * c - d + 34) % 7;
    int j = d + e + 114;
    int month = j / 31;
    int day = j % 31 + 1;
    if (gregorian_dates)
        return jdn_to_gregorian(julian_date_to_jdn(year, month, day));
    return {year, month, day};
}

double prob_julian_easter(int day) {
    if (day < 0 || day >= 43) return 0.0;
    return static_cast<double>(JD[day]) / 1064.0;  // TWOCYCLE = CYCLE << 1
}

double prob_easter(int day) {
    if (day < 0 || day >= 35) return 0.0;
    if (day < 6) return static_cast<double>(day + 1) / 206.71413;
    if (day < 28) return 0.03386319067786996;
    return (static_cast<double>(35 - day) + 0.53059) / 206.71413;
}

// ---------------------------------------------------------------- family --  //

std::vector<std::vector<int>> td_count(const calendar_domain& domain) {
    validate_domain(domain);
    require_trading_frequency(domain, "trading-day generation");
    const int r = ratio_of(domain.frequency);
    std::vector<std::vector<int>> out(7, std::vector<int>(
                                           static_cast<std::size_t>(domain.length), 0));
    for (int i = 0; i < domain.length; ++i) {
        int m0 = domain.month + i * r;
        int y = domain.year + (m0 - 1) / 12;
        int m = (m0 - 1) % 12 + 1;
        for (int mm = m; mm < m + r; ++mm) {
            int yy = y + (mm - 1) / 12;
            int mo = (mm - 1) % 12 + 1;
            for (int d = 1; d <= month_len(yy, mo); ++d)
                ++out[weekday_mon(yy, mo, d)][static_cast<std::size_t>(i)];
        }
    }
    return out;
}

std::vector<double> mean_days_yearly(int frequency) {
    if (frequency != 12)
        throw std::invalid_argument(
            "calendar_generation: mean_days_yearly supports frequency 12 only");
    std::vector<double> out(12);
    for (int i = 0; i < 12; ++i) out[static_cast<std::size_t>(i)] = MDAYS[i] / 7.0;
    return out;
}

std::vector<double> mean_days(const calendar_domain& domain) {
    validate_domain(domain);
    const int r = ratio_of(domain.frequency);
    std::vector<double> out(static_cast<std::size_t>(domain.length));
    for (int i = 0; i < domain.length; ++i) {
        int m0 = domain.month + i * r;
        int m = (m0 - 1) % 12 + 1;
        double acc = 0.0;
        for (int mm = m; mm < m + r; ++mm) acc += MDAYS[((mm - 1) % 12)];
        out[static_cast<std::size_t>(i)] = acc / 7.0;
    }
    return out;
}

std::vector<std::vector<double>> trading_days_contrast(day_clustering clustering,
                                                       const calendar_domain& domain) {
    const cluster_info& c = cluster_at(clustering);
    std::vector<std::vector<int>> counts = td_count(domain);  // gates frequencies
    std::vector<std::vector<int>> members(static_cast<std::size_t>(c.groups));
    for (int w = 0; w < 7; ++w)
        members[static_cast<std::size_t>(c.def[w])].push_back(w);
    const int ref = 0;
    const double xd = static_cast<double>(members[static_cast<std::size_t>(ref)].size());
    std::vector<std::vector<double>> cols;
    for (int g = 1; g < c.groups; ++g) {
        const auto& mg = members[static_cast<std::size_t>(g)];
        const double size = static_cast<double>(mg.size());
        std::vector<double> col(static_cast<std::size_t>(domain.length));
        for (int i = 0; i < domain.length; ++i) {
            int sum_ref = 0;
            for (int w : members[static_cast<std::size_t>(ref)])
                sum_ref += counts[static_cast<std::size_t>(w)][static_cast<std::size_t>(i)];
            int sum_g = 0;
            for (int w : mg) sum_g += counts[static_cast<std::size_t>(w)][static_cast<std::size_t>(i)];
            const double x = static_cast<double>(sum_ref) / xd;  // divide FIRST
            col[static_cast<std::size_t>(i)] =
                static_cast<double>(sum_g) - size * x;
        }
        cols.push_back(std::move(col));
    }
    return cols;
}

std::vector<double> length_of_period_column(calendar_lop type,
                                            const calendar_domain& domain) {
    validate_domain(domain);
    const std::size_t n = static_cast<std::size_t>(domain.length);
    std::vector<double> col(n, 0.0);
    switch (type) {
        case calendar_lop::none:
            return col;
        case calendar_lop::leap_year: {
            // LPFactory.lp: February (index k = 1 for monthly) or the first
            // period of each year (k = 0 otherwise) strikes every period;
            // the first leap year and every 4 periods from it strike 0.75
            // (4-year-assumption quirk; bounds-free search keeps going past
            // the range, mirroring the arithmetic TsDomain.get).
            const int f = domain.frequency;
            const int k = f == 12 ? 1 : 0;
            int m0 = (domain.month - 1) % f;
            int p = k - m0;
            if (p < 0) p += f;
            for (int i = p; i < domain.length; i += f) col[static_cast<std::size_t>(i)] = -0.25;
            long idx = static_cast<long>(domain.month - 1) + static_cast<long>(p) * ratio_of(f);
            int q = p;
            int y = domain.year + static_cast<int>(idx / 12);
            while (!is_leap(y)) {
                q += f;
                ++y;
            }
            for (int i = q; i < domain.length; i += 4 * f) col[static_cast<std::size_t>(i)] = 0.75;
            return col;
        }
        case calendar_lop::length_of_period: {
            const double c = 365.25 / static_cast<double>(domain.frequency);
            const int r = ratio_of(domain.frequency);
            for (int i = 0; i < domain.length; ++i) {
                int m0 = domain.month + i * r;
                int y = domain.year + (m0 - 1) / 12;
                int m = (m0 - 1) % 12 + 1;
                int days = 0;
                for (int mm = m; mm < m + r; ++mm)
                    days += month_len(y + (mm - 1) / 12, (mm - 1) % 12 + 1);
                col[static_cast<std::size_t>(i)] = static_cast<double>(days) - c;
            }
            return col;
        }
    }
    throw std::invalid_argument("calendar_generation: unknown length-of-period type");
}

std::vector<double> easter_column(int duration, easter_correction correction,
                                  int end_position, const calendar_domain& domain) {
    validate_domain(domain);
    require_trading_frequency(domain, "easter generation");
    if (duration < 1 || duration > 28)
        throw std::invalid_argument(
            "calendar_generation: easter duration must be in [1, 28]");
    if (correction == easter_correction::pre_computed && duration > 25)
        throw std::invalid_argument(
            "calendar_generation: PreComputed correction requires a duration in [1, 25]");
    const int f = domain.frequency;
    const int n = domain.length;
    const std::size_t un = static_cast<std::size_t>(n);
    std::vector<double> col(un, 0.0);
    const int k12 = 12 / f;
    int year = domain.year;
    int month = domain.month;
    if (month > 4) {
        month -= 12;
        ++year;
    }
    int pos = (4 - month) / k12;
    if (pos > n) return col;
    const double durd = static_cast<double>(duration);
    while (pos <= n) {
        calendar_date e = easter_date(year);
        int day = e.day;
        int mon = e.month;
        if (end_position == -1) {
            --day;
        } else if (end_position == 1) {
            if (day == 31) {
                day = 1;
                ++mon;
            } else {
                ++day;
            }
        }
        double c1 = 0.0, c2 = 0.0;
        if (correction == easter_correction::simple) {
            if (mon == 3 || day == 0) {
                c1 = 0.5;
                c2 = -0.5;
            } else if (day >= duration) {
                c1 = -0.5;
                c2 = 0.5;
            } else {
                c1 = (durd - day) / durd - 0.5;
                c2 = -c1;
            }
        } else {
            double m25 = 0.0, m27 = 0.0;
            if (correction == easter_correction::pre_computed) {
                m25 = MEANS_MAR[25 - duration];
                m27 = MEANS_APR[25 - duration];
            } else if (correction == easter_correction::theoretical) {
                m25 = (15.0 + durd) / 59.06118;
                m27 = 1 - m25;
            }
            if (mon == 3) {
                c1 = 1 - m25;
                c2 = -m27;
            } else if (day >= duration) {
                c1 = -m25;
                c2 = 1 - m27;
            } else {
                c1 = (durd - day) / durd - m25;
                c2 = day / durd - m27;
            }
        }
        if (pos - 1 >= 0) col[static_cast<std::size_t>(pos - 1)] = c1;
        if (pos < n) col[static_cast<std::size_t>(pos)] = c2;
        ++year;
        pos += f;
    }
    return col;
}

std::vector<double> julian_easter_column(int duration, bool gregorian_dates,
                                         const calendar_domain& domain) {
    validate_domain(domain);
    require_trading_frequency(domain, "julian easter generation");
    if (duration < 1 || duration > 28)
        throw std::invalid_argument(
            "calendar_generation: julian easter duration must be in [1, 28]");
    const int f = domain.frequency;
    const int n = domain.length;
    std::vector<double> col(static_cast<std::size_t>(n), 0.0);
    const double tot = static_cast<double>(1064 * duration);
    const double durd = static_cast<double>(duration);
    const long dstart = days_from_civil(domain.year, static_cast<unsigned>(domain.month), 1);
    int em0 = domain.month + n * ratio_of(f);
    const long dend =
        days_from_civil(domain.year + (em0 - 1) / 12,
                        static_cast<unsigned>((em0 - 1) % 12 + 1), 1);
    for (int year = domain.year;
         days_from_civil(year, 3, 1) <= dend; ++year) {
        calendar_date ec = julian_easter_date(year, gregorian_dates);
        const long edays = days_from_civil(ec.year, static_cast<unsigned>(ec.month),
                                           static_cast<unsigned>(ec.day));
        if (!(dstart < edays)) continue;
        const long w0 = edays - duration;  // window [e - dur, e - 1]
        int x19 = static_cast<int>(days_from_civil(year, 4, 1) - w0);
        if (x19 < 0) x19 = 0;              // mar.end() is EXCLUSIVE (April 1)
        int x20 = static_cast<int>(edays - days_from_civil(year, 5, 1));
        if (x20 < 0) x20 = 0;
        else if (x20 > duration) x20 = duration;
        const int x21 = duration - x19 - x20;
        const int months[3] = {3, 4, 5};
        const int x[3] = {x19, x21, x20};
        const int* tables[3] = {C_MAR2, C_APR2, C_MAY2};
        for (int t = 0; t < 3; ++t) {
            const int idx = pos_of(domain, year, months[t]);
            if (idx >= 0 && idx < n) {
                col[static_cast<std::size_t>(idx)] +=
                    static_cast<double>(x[t]) / durd -
                    static_cast<double>(tables[t][duration - 1]) / tot;
            }
        }
    }
    return col;
}

std::vector<double> ramp_column(int start_year, int start_month, int end_year,
                                int end_month, const calendar_domain& domain) {
    validate_domain(domain);
    for (int m : {start_month, end_month})
        if (m < 1 || m > 12)
            throw std::invalid_argument("calendar_generation: ramp month out of range");
    for (int y : {start_year, end_year})
        if (y < 1 || y > 9999)
            throw std::invalid_argument(
                "calendar_generation: ramp year out of range [1, 9999]");
    const int n = domain.length;
    std::vector<double> col(static_cast<std::size_t>(n), 0.0);
    // TsPeriod.of(unit, date) truncates to the period containing the month.
    const int i1 = pos_of(domain, start_year, start_month);
    const int i2 = pos_of(domain, end_year, end_month);
    if (i2 <= i1) return col;  // reversed/equal ramp: all zeros
    if (i1 >= 0 && i1 >= n - 1) {
        std::fill(col.begin(), col.end(), -1.0);
        return col;
    }
    if (i1 >= 0)
        for (int i = 0; i <= i1; ++i) col[static_cast<std::size_t>(i)] = -1.0;
    if (i2 < n) {
        if (i2 <= 0) {
            std::fill(col.begin(), col.end(), 0.0);
            return col;
        }
        for (int i = i2; i < n; ++i) col[static_cast<std::size_t>(i)] = 0.0;
    }
    const int lo = std::max(i1 + 1, 0);
    const int hi = std::min(i2, n);
    const double den = static_cast<double>(i2 - i1);
    for (int i = lo; i < hi; ++i)
        col[static_cast<std::size_t>(i)] =
            static_cast<double>(i - i1) / den - 1.0;
    return col;
}

std::vector<double> user_column(const std::vector<double>& data, int data_year,
                                int data_month, const calendar_domain& domain) {
    validate_domain(domain);
    const int f = domain.frequency;
    if (f == 12 && (data_month < 1 || data_month > 12))
        throw std::invalid_argument(
            "calendar_generation: user data month not aligned with the domain frequency");
    if (f == 4 && data_month != 1 && data_month != 4 && data_month != 7 &&
        data_month != 10)
        throw std::invalid_argument(
            "calendar_generation: user data month not aligned with the domain frequency");
    if (f == 1 && data_month != 1)
        throw std::invalid_argument(
            "calendar_generation: user data month not aligned with the domain frequency");
    if (data_year < 1 || data_year > 9999)
        throw std::invalid_argument(
            "calendar_generation: user data year out of range [1, 9999]");
    if (data.empty())
        throw std::invalid_argument("calendar_generation: user data must be non-empty");
    const int n = domain.length;
    const int dlen = static_cast<int>(data.size());
    // TsVariablesFactory: pos = dataStart.until(p) (negative when the data
    // starts at/after the domain), end = pos + n; overlap [from, to) receives
    // data[pos, end); everything else keeps the zero fill.
    const int dpos = pos_of(domain, data_year, data_month);
    int pos = -dpos;
    int end = pos + n;
    int from = 0, to = n;
    if (pos < 0) {
        from = -pos;
        pos = 0;
    }
    if (end > dlen) {
        to -= end - dlen;
        end = dlen;
    }
    std::vector<double> col(static_cast<std::size_t>(n), 0.0);
    if (from < to)
        std::copy(data.begin() + pos, data.begin() + end, col.begin() + from);
    return col;
}

// ---------------------------------------------------------------- names ----

const char* trading_days_tag(day_clustering clustering) {
    return cluster_at(clustering).tag;
}

const char* trading_days_name(day_clustering clustering) {
    return cluster_at(clustering).name;
}

const char* trading_days_default_name(int coefficient_count) {
    if (coefficient_count < 1 || coefficient_count > 8)
        throw std::invalid_argument(
            "calendar_generation: defaultName supports 1..8 (captured range)");
    return "td";
}

const char* trading_days_variable_name() { return "td"; }
const char* length_of_period_variable_name() { return "lp"; }

std::string easter_description(int duration) {
    if (duration < 1 || duration > 28)
        throw std::invalid_argument(
            "calendar_generation: easter duration must be in [1, 28]");
    return "easter [" + std::to_string(duration) + "]";
}

std::string julian_easter_description(int duration) {
    if (duration < 1 || duration > 28)
        throw std::invalid_argument(
            "calendar_generation: julian easter duration must be in [1, 28]");
    return "julian easter [" + std::to_string(duration) + "]";
}

const char* length_of_period_description(calendar_lop type) {
    switch (type) {
        case calendar_lop::none:
        case calendar_lop::length_of_period:
            return "Length of period";
        case calendar_lop::leap_year:
            return "Leap year";
    }
    throw std::invalid_argument("calendar_generation: unknown length-of-period type");
}

std::string ramp_description(int start_year, int start_month, int end_year,
                             int end_month, int frequency) {
    if (frequency != 12)
        throw std::invalid_argument(
            "calendar_generation: ramp_description supports monthly domains only");
    for (int m : {start_month, end_month})
        if (m < 1 || m > 12)
            throw std::invalid_argument("calendar_generation: ramp month out of range");
    for (int y : {start_year, end_year})
        if (y < 1 || y > 9999)
            throw std::invalid_argument(
                "calendar_generation: ramp year out of range [1, 9999]");
    char buf[64];
    std::snprintf(buf, sizeof(buf), "ramp %04d-%02d / %04d-%02d", start_year,
                  start_month, end_year, end_month);
    return std::string(buf);
}

const std::vector<trading_days_type_info>& trading_days_types() {
    static const std::vector<trading_days_type_info> T = {
        {"NONE", 0, "null", "null", "null"},
        {"TD7", 7, "[mon, tue, wed, thu, fri, sat, sun]",
         "[mon, tue, wed, thu, fri, sat]", "Mon/Tue/Wed/Thu/Fri/Sat/Sun"},
        {"TD4", 4, "[mon-thu, fri, sat, sun]", "[mon-thu, fri, sat]",
         "Mon-Thu/Fri/Sat/Sun"},
        {"TD3c", 3, "[mon-thu, fri-sat, sun]", "[mon-thu, fri-sat]",
         "Mon-Thu/Fri-Sat/Sun"},
        {"TD3", 3, "[mon-fri, sat, sun]", "[mon-fri, sat]", "Mon-Fri/Sat/Sun"},
        {"TD2", 2, "[week days, week-end]", "[week days]",
         "Week days/Week-end"},
        {"TD2c", 2, "[mon-sat, sun]", "[mon-sat]", "Mon-Sat/Sun"},
        {"TD2d", 2, "null", "[mon-thu]", "null"},
        {"TD4c", 4, "[mon, tue-fri, sat, sun]", "[mon, tue-fri, sat]",
         "Mon/Tue-Fri/Sat/Sun"},
        {"TD6", 6, "[mon, tue, wed, thu, fri, week-end]",
         "[mon, tue, wed, thu, fri]", "Mon/Tue/Wed/Thu/Fri/Week-end"},
        {"TDuser", -1, "null", "null", "null"},
    };
    return T;
}

}  // namespace kalos
