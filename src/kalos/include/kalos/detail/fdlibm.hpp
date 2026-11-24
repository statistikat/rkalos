#pragma once

#include <stdexcept>

namespace kalos::detail {

// Transcription of the pinned JDK 21 `java.lang.StrictMath`/`Math` transcendental
// kernels (java.lang.FdLibm.Exp.compute, java.lang.FdLibm.Log.compute) as dumped
// from the pinned runtime's bytecode. The jdplus machine calls Math.exp and
// Math.log, which delegate to these kernels, so the qualified numerics reproduce
// them operation for operation: the C library's exp/log differ by an ulp at
// scattered arguments (for example the scaling factor exp(ldet/n) of the
// deviances machine), which is enough to steer the estimation chain into a
// different basin. The pinned capture JVM additionally replaces Math.exp/log
// with build-specific StubRoutines::dexp/dlog assembly at hot call sites; the
// corpus is therefore captured with those intrinsics disabled (see the item
// 21-D contract) so that Math.exp/log equal the FdLibm bytecode semantics.
double fdlibm_exp(double x);
double fdlibm_log(double x);
// StrictMath.pow(x, 1.0/3.0) for the portable residual normality profile.
double fdlibm_power_one_third(double x);

// CI-FDLIBM routing kernels (provenance table, pinned runtime
// OpenJDK 21.0.11+10-1-24.04.2-Ubuntu, bytecode dumps under
// tmp_opencode/ci-fdlibm/dumps/, golden-proof records under
// tmp_opencode/ci-fdlibm/goldens/):
//
//   kernel         target of record                                  JVM relation
//   fdlibm_expm1   java.lang.FdLibm$Expm1.compute                    StrictMath.expm1 dispatch target
//                                                                                   (Math.expm1 -> StrictMath.expm1; no
//                                                                                   HotSpot intrinsic exists); bit-exact
//                                                                                   vs double-captured JVM goldens on
//                                                                                   4.4M grid + all corpus-reached args
//   fdlibm_tanh    java.lang.FdLibm$Tanh.compute                     Math.tanh -> StrictMath.tanh ->
//                                                                                   this kernel (JDK 9+ rewrote the
//                                                                                   hyperbolic kernel: it is expm1-based,
//                                                                                   NOT fdlibm-5.x e_tanh.c); bit-exact
//                                                                                   on 4.5M grid + corpus args
//   fdlibm_log1p   java.lang.FdLibm$Log1p.compute                    StrictMath.log1p dispatch target;
//                                                                                   bit-exact on 8.5M grid + corpus args
//   fdlibm_atanh   fdlibm-5.x __ieee754_atanh (netlib e_atanh.c      NO JVM counterpart exists: neither
//                  1.3 95/01/18) + fdlibm_log1p                      java.lang.Math nor StrictMath has
//                                                                                   atanh in the pinned JDK. The kernel is
//                                                                                   the fdlibm-5.x sibling of the vendored
//                                                                                   kernels, chosen for machine-
//                                                                                   independence of the reflection-
//                                                                                   coordinate map; audited <= 1 ulp vs
//                                                                                   glibc (2M grid + corpus args),
//                                                                                   specials + tanh round-trip exact.
//   fdlibm_cos     java.lang.FdLibm$Cos.compute (+$Sin, $RemPio2)    Bytecode (StrictMath) semantics; bit-
//   fdlibm_sin     java.lang.FdLibm$Sin.compute                      exact on 10M-point grids (incl. |x|
//                                                                                   up to 360 through the staged Cody-Waite
//                                                                                   fix-up branches) + all corpus args.
//
// Boundary: the mid-range Cody-Waite path (RemPio2, ix <= 0x413921FB, i.e.
// |x| <= ~1.646e6) is complete; the huge-argument Payne-Hanek path
// ($KernelRemPio2) is NOT transcribed and throws std::domain_error.
//
// IMPORTANT — Math.cos/sin on this runtime (probe 2026-11-21, supersedes an
// earlier wrong note here): the _dcos/_dsin C2 stubs ARE live on this Intel
// CPU (HotSpot 21 x86 stubs are generated code, not XOP). Their bits differ
// from the FdLibm bytecode by 1 ulp at scattered arguments (337/14753 cos,
// 399/9507 sin of the corpus-reached set), and they agree with glibc instead
// (19/16 disagreements) — HotSpot's libm stubs share the IBM-Accurate-
// Portable-Math-Library lineage that glibc and the macOS libm also derive
// from, which is why raw std::cos/std::sin reproduce the captured corpora on
// every qualified platform. CONSEQUENCE: the routed surface deliberately
// EXCLUDES cos/sin (staying on host libm matches the corpus bit lineage; the
// bytecode kernels above are vendored for their other JVM-fidelity roles and
// the golden-proof harness). For golden captures of any Math.* semantics use
// StrictMath (or -XX:DisableIntrinsic=_dexp -XX:DisableIntrinsic=_dlog
// -XX:DisableIntrinsic=_dcos -XX:DisableIntrinsic=_dsin): expm1/log1p/tanh
// have no intrinsic and Math==StrictMath for them bit-for-bit.
double fdlibm_expm1(double x);
double fdlibm_tanh(double x);
double fdlibm_log1p(double x);
double fdlibm_atanh(double x);  // domain: |x| <= 1, see table above
double fdlibm_cos(double x);    // throws std::domain_error when hi_word(|x|) > 0x413921FB (~1.646e6)
double fdlibm_sin(double x);    // throws std::domain_error when hi_word(|x|) > 0x413921FB (~1.646e6)

} // namespace kalos::detail
