#pragma once

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

} // namespace kalos::detail
