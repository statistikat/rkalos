#pragma once

// Host-lineage transcendental kernels: a vendored, machine-independent
// transcription of the exact glibc-2.39 dbl-64 kernels that the qualification
// host's libm.so.6 executes. These exist so that every corpus-reachable
// bit-claim site computes the SAME result bits on any little-endian IEEE-754
// target, independent of the C library (CI-FDLIBM remediation, LEG A of
// tmp_opencode/kalos-ci3-fix/fix-proposal.PROPOSED.md; ledger:
// docs/contracts/x13-automatic-pipeline-contract.md, CI-FDLIBM section).
//
// Unlike the fdlibm_* kernels in fdlibm.hpp — vendored to reproduce the
// pinned JVM's bytecode kernels — these reproduce the HOST LIBM lineage that
// the captured corpora were qualified against. That lineage is the IBM
// Accurate Portable Math Library (sin/cos/atan2/acos) plus the ARM
// table-driven kernels (exp/log/pow) and the Borges MyHypot3 kernel (hypot),
// as shipped in Ubuntu 24.04's glibc 2.39-0ubuntu8.8. Empirically this is
// also the lineage HotSpot's _dcos/_dsin stubs and the macOS libm derive
// from (fdlibm.hpp probe note, 2026-11-21), which is why raw std::cos/sin
// reproduced the corpora on every qualified platform pre-vendoring.
//
// Provenance (retrieved 2026-11-22 from
// https://raw.githubusercontent.com/bminor/glibc/glibc-2.39/sysdeps/ieee754/dbl-64/,
// official GNU mirror; glibc tag glibc-2.39):
//
//   file                   role                                            SHA256 (file text)
//   s_sin.c                __sin/__cos + do_sin/do_cos/reduce_sincos       d607bcb0823bb793ecf917546266559aa6f3da6504719ecb90db8b6cb4bb8b1b
//   branred.c              __branred (huge-argument reduction)             e697adc29fed48602ffd3c790f90bda8da1db7bc5e135d39f9880a0f7480d6ad
//   branred.h              branred constants + toverp[75]                  6b67701f64741b716ba8915569dcb1406095e1bf60596a841df2985464b3c3f2
//   usncs.h                sin/cos polynomial + reduction constants        b69a58a10c9f2f2846eb16bd1c4f9fac8f4b2f4e33873348d1b5cb3f1a9c2c72
//   sincostab.c            __sincostab table (via generator)               3233d23c6085ac87d47ce1fd7471d88533213f6271ba333e664df30f3172ab85
//   e_exp.c                __exp                                           68867b659baecfaa37aa81930d74fd91dc44333ed979ac630a1c06dc5b9d4ce8
//   e_exp_data.c           __exp_data table (via generator)                9dee45e089bc8fac36a79da80889c7d784793aae7a00743a32066bd49dbf3434
//   e_log.c                __log                                           a645fd8877ca71441cbf4268ab1bdd2ca38f3054ebfe15d75e820742f75770ae
//   e_log_data.c           __log_data table (via generator)                614b871035d0511f2ebfdbc5f7305f1974d2087ae3fba6aeed0f42dd5f091ba5
//   e_pow.c                __pow (own log_inline/exp_inline variants)      68c6da4de5b99d1c0500d6b5078a702274d22cd8fb77a47a6a21ca025b6455e2
//   e_pow_log_data.c       __pow_log_data table (via generator)            c64bee05138417398b2d2b2a08c67d96cf22a135c1b6a72509cdadb4edf633be
//   e_hypot.c              __hypot (Borges MyHypot3)                       5bb30147cb13eade6e28eddac0a69b04de71d2850090166e36c94236ac9c514e
//   e_atan2.c              __ieee754_atan2 (IBM APML uatan2)               d0dcea030cfa9054d5daeecd19531bc18ee33fa93714b28bae4dca474d290463
//   atnat2.h               atan2 constants                                 7a6f80ed62403a984ea78592b6e98b39361999855971198cd0f454294d1d9db9
//   uatan.tbl              cij[241][7] table (via generator)               8072e14d43f1b897ab7013b70954e18b113ef2fcac942e8a263a579d4baf531a
//   e_asin.c               __ieee754_asin/__ieee754_acos (IBM APML)        cf56fcbf627fd8d57ca280e95ec891b978c70ae54f0f5929cc017003a3256eeb
//   uasncs.h               asin/acos constants                             e6dc99829e93c7208103ea62bbbd8a211560afa8107509631b1ba9690635b56e
//   asincos.tbl            asncs table (via generator)                     067c3b7047ab28b087262ba26a839ad23bf39236a843897306d5c828d69904f0
//   root.tbl               inroot table (via generator)                    f3c9d14fdf380e8db2098a183266721cd800ddeb4745ef14f2a491c35b17f5eb
//   powtwo.tbl             powtwo table (via generator)                    e66e94adcd700010aa6ea5407130486ef4dba610311c1aab74bda7952bd67f53
//   dla.h                  EADD/ESUB/EMULV double-length macros            466e4b6c0092f6d104e310bb4f3e79e97deaeccccd94e5a21f6dc7b4d71367c4
//   mydefs.h               mynumber/int4/HIGH_HALF plumbing                84a11183a0d202df5ccdbe59a53eaec131c28d36e53db29bd9118a1a65dbf9bf
//   math_config.h          shared plumbing + table struct layouts          90c710c89759618256073f9dbef7a5cc662f4a1cfe94df1ea3d7b184ae96845c
//
// The table DATA lives in src/hostlineage_tables.cpp, generated mechanically
// from the pinned *.tbl/sincostab.c/e_*_data.c inputs by
// tools/gen_hostlineage_tables.py (digests pinned in
// tools/hostlineage-glibc-sha256.json; the generator refuses unexpected
// guards/counts). tools/verify_hostlineage.cpp re-checks every transcribed
// table element against values parsed independently from the original files
// AND proves bit-equality of every kernel below against the host libm over
// dense argument grids and the FULL corpus-reached argument census recorded
// under the shim (tmp_opencode/kalos-ci4-fix/rec2/).
//
// Lineage decisions (each validated by the proof tool on the qualification
// host, which IS the reference: glibc 2.39, x86-64 baseline):
//   * TOINT_INTRINSICS = 0 (Shift-path argument reduction in exp/pow-exp).
//     glibc never sets this on x86-64; the roundtoint branch is not
//     transcribed.
//   * __FP_FAST_FMA not defined (non-FMA branches everywhere: exp/log/pow
//     table selection, hypot MyHypot3 correction, dla.h EMULV). Ubuntu's
//     libm.so.6 contains NO STT_GNU_IFUNC symbols (nm -D audit 2026-11-22),
//     i.e. one non-FMA variant of each kernel is what any process on the
//     qualified hosts executes; the FMA branches are #if'd out.
//   * WANT_ROUNDING = 1, WANT_ERRNO = 0. Return-value BITS are pinned; the
//     errno/FE_* exception side effects of the glibc kernels are
//     deliberately DROPPED (kalos never reads them; math_force_eval /
//     __set_errno sites are omitted where they cannot affect returned bits).
//   * SET_RESTORE_ROUND(_53BIT) is a no-op: the kernels assume
//     round-to-nearest, which kalos never changes and the whole numeric
//     surface already assumes.
//   * math_narrow_eval is the identity, valid because FLT_EVAL_METHOD == 0
//     is static_asserted on the target.
//   * __math_uflow/oflow/divzero/invalid are the bit-equivalent constants of
//     glibc's s_math_*.c (±0, ±inf, ±inf, (x-x)/(x-x)); check_oflow/
//     check_uflow are pass-throughs (WANT_ERRNO = 0).
//   * Statement-expression macros (TAYLOR_SIN, SINCOS_TABLE_LOOKUP, EADD,
//     ESUB, EMULV) became inline reference-output functions with verbatim
//     bodies (-Wpedantic cleanliness); bodies/operation order are
//     character-faithful.
//   * Only the LITTLE_ENDI table branches are vendored (the proof oracle is
//     the little-endian host; big-endian targets fail closed via
//     static_assert in src/hostlineage_internal.hpp).
//
// Licensing FLAG (open maintainer issue, recorded 2026-11-22): this repo is
// MIT, but the vendored files are glibc LGPL-2.1-or-later (IBM APML files:
// LGPL-2.1+; ARM table kernels: LGPL-2.1+; Borges hypot: LGPL-3.0+ with
// linking exception per upstream). Vendoring as static-source-with-provenance
// is the CI-remediation scope agreed for this slice; a licensing sign-off or
// a dynamic-provider alternative is a follow-up for the maintainers.

#include <cmath>

namespace kalos::detail {

// Routed bit-claim surfaces (CI-FDLIBM ledger, LEG A flips):
//   sin/cos  - ucarima seasonal harmonics (ucarima_engine.hpp,
//              seats_ucarima_model.hpp, seats_linear_filters.hpp)
//   exp/log  - x13_automatic/x13_fit scale exp, x13_allocation log-additive
//              projection, x11_log_additive log
//   pow      - preparation/outlier_scoring/ucarima-engine integer powers
//   hypot    - diffuse_state_space, regression_likelihood,
//              finite_state_factor, sparse_history, ssq_levmar (bit-exact
//              claim site), and the seats linear-filter path
//   atan2    - ucarima_engine.hpp:90 decomposition phase
//   acos     - ucarima_engine.hpp:1143 seasonal period recovery
double glibc_sin(double x);
double glibc_cos(double x);
// x86-64 HotSpot cosine profile for bounded SEATS spectrum frequencies.
double hotspot_spectrum_cos(double x);
double glibc_exp(double x);
double glibc_log(double x);
double glibc_pow(double x, double y);
double glibc_hypot(double x, double y);
double glibc_atan2(double y, double x);
double glibc_acos(double x);

// Exposed for proof completeness only (same transcription unit as acos;
// nothing in kalos routes to it — the ledger has no asin sites).
double glibc_asin(double x);

}  // namespace kalos::detail
