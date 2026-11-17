#include "kalos/detail/lapack.hpp"
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <type_traits>
static_assert(sizeof(int) == 4 && std::is_same<La_INT, int>::value &&
              std::is_same<BLAS_INT, int>::value, "Kalos requires R's 32-bit BLAS/LAPACK ABI");
namespace kalos { namespace lapack {
void dgelsy(const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb, int* jpvt, const double* rcond, int* rank, double* work, const int* lwork, int* info) {
    F77_CALL(dgelsy)(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
}
void dgeev(const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info) {
    F77_CALL(dgeev)(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info FCONE FCONE);
}
void dgeqrf(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info) {
    F77_CALL(dgeqrf)(m, n, a, lda, tau, work, lwork, info);
}
void dgesvd(const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info) {
    F77_CALL(dgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info FCONE FCONE);
}
void dpocon(const char* uplo, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info) {
    F77_CALL(dpocon)(uplo, n, a, lda, anorm, rcond, work, iwork, info FCONE);
}
void dpotrf(const char* uplo, const int* n, double* a, const int* lda, int* info) {
    F77_CALL(dpotrf)(uplo, n, a, lda, info FCONE);
}
void dpotrs(const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda, double* b, const int* ldb, int* info) {
    F77_CALL(dpotrs)(uplo, n, nrhs, a, lda, b, ldb, info FCONE);
}
void dsyev(const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* w, double* work, const int* lwork, int* info) {
    F77_CALL(dsyev)(jobz, uplo, n, a, lda, w, work, lwork, info FCONE FCONE);
}
void dlarfg(const int* n, double* alpha, double* x, const int* incx, double* tau) {
    F77_CALL(dlarfg)(n, alpha, x, incx, tau);
}
void dgemm(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, const double* beta, double* c, const int* ldc) {
    F77_CALL(dgemm)(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc FCONE FCONE);
}
const char* provider_name() noexcept { return "R-configured BLAS/LAPACK"; }
const char* provider_version() noexcept { return "configured by R"; }
}}
