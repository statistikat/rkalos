// Direct LAPACK singular-value probes for the identification-fixture tolerance
// study (see tools/svd-tolerance-evidence.md). They call the same F77_CALL
// boundary as src/lapack_r.cpp so that R's svd() target, DGESVD and DGESDD can be
// attributed to the same shared library. Not part of the package build.
//
//   PKG_LIBS="$(R CMD config LAPACK_LIBS) $(R CMD config BLAS_LIBS) $(R CMD config FLIBS)" \
//     R CMD SHLIB -o svd-tolerance-probe.so tools/svd-tolerance-probe.cpp
#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include <vector>

// Singular values only, exactly as kalos::singular_values() requests them.
extern "C" SEXP probe_dgesvd(SEXP matrix) {
    if (TYPEOF(matrix) != REALSXP || !Rf_getAttrib(matrix, R_DimSymbol))
        error("matrix argument must be a numeric matrix");
    const int m = INTEGER(getAttrib(matrix, R_DimSymbol))[0];
    const int n = INTEGER(getAttrib(matrix, R_DimSymbol))[1];
    const int small = m < n ? m : n;
    SEXP copy = PROTECT(allocVector(REALSXP, (R_xlen_t) m * n));
    memcpy(REAL(copy), REAL(matrix), sizeof(double) * (size_t) (m * n));
    double* a = REAL(copy);
    std::vector<double> values(small);
    double query = 0;
    int one = 1, info = 0, lwork = -1;
    F77_CALL(dgesvd)("N", "N", &m, &n, a, &m, values.data(), a, &one, a, &one, &query, &lwork, &info FCONE FCONE);
    lwork = (int) query;
    std::vector<double> work(lwork);
    F77_CALL(dgesvd)("N", "N", &m, &n, a, &m, values.data(), a, &one, a, &one, work.data(), &lwork, &info FCONE FCONE);
    SEXP out = PROTECT(allocVector(REALSXP, small));
    for (int i = 0; i < small; ++i) REAL(out)[i] = values[i];
    setAttrib(out, install("info"), ScalarInteger(info));
    UNPROTECT(2);
    return out;
}

// Divide-and-conquer driver, singular values only.
extern "C" SEXP probe_dgesdd(SEXP matrix) {
    if (TYPEOF(matrix) != REALSXP || !Rf_getAttrib(matrix, R_DimSymbol))
        error("matrix argument must be a numeric matrix");
    const int m = INTEGER(getAttrib(matrix, R_DimSymbol))[0];
    const int n = INTEGER(getAttrib(matrix, R_DimSymbol))[1];
    const int small = m < n ? m : n;
    SEXP copy = PROTECT(allocVector(REALSXP, (R_xlen_t) m * n));
    memcpy(REAL(copy), REAL(matrix), sizeof(double) * (size_t) (m * n));
    double* a = REAL(copy);
    std::vector<double> values(small);
    std::vector<int> iwork(8 * small);
    double query = 0;
    int one = 1, info = 0, lwork = -1;
    F77_CALL(dgesdd)("N", &m, &n, a, &m, values.data(), a, &one, a, &one, &query, &lwork, iwork.data(), &info FCONE);
    lwork = (int) query;
    std::vector<double> work(lwork);
    F77_CALL(dgesdd)("N", &m, &n, a, &m, values.data(), a, &one, a, &one, work.data(), &lwork, iwork.data(), &info FCONE);
    SEXP out = PROTECT(allocVector(REALSXP, small));
    for (int i = 0; i < small; ++i) REAL(out)[i] = values[i];
    setAttrib(out, install("info"), ScalarInteger(info));
    UNPROTECT(2);
    return out;
}

extern "C" void R_init_svd_tolerance_probe(DllInfo* dll) {
    static const R_CallMethodDef calls[] = {
        {"probe_dgesvd", (DL_FUNC) &probe_dgesvd, 1},
        {"probe_dgesdd", (DL_FUNC) &probe_dgesdd, 1},
        {NULL, NULL, 0}
    };
    R_registerRoutines(dll, NULL, calls, NULL, NULL);
}
