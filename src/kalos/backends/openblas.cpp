#include "kalos/detail/lapack.hpp"
#include <cblas.h>
#include <lapacke.h>
#include <stdexcept>
#include <type_traits>
static_assert(sizeof(int) == 4 && std::is_same<lapack_int, int>::value,
              "Kalos requires the 32-bit LP64 OpenBLAS integer ABI");
namespace kalos { namespace lapack {
namespace {
int fortran_info(int info) { return info < 0 ? info + 1 : info; }
CBLAS_TRANSPOSE transpose(char c) {
    if (c == 'N') return CblasNoTrans;
    if (c == 'T') return CblasTrans;
    if (c == 'C') return CblasConjTrans;
    throw std::invalid_argument("Unsupported BLAS transpose flag");
}
}
void dgelsy(const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb, int* jpvt, const double* rcond, int* rank, double* work, const int* lwork, int* info) {
    *info = fortran_info(LAPACKE_dgelsy_work(LAPACK_COL_MAJOR, *m, *n, *nrhs, a, *lda, b, *ldb, jpvt, *rcond, rank, work, *lwork));
}
void dgeqrf(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info) {
    *info = fortran_info(LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, *m, *n, a, *lda, tau, work, *lwork));
}
void dgesvd(const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info) {
    *info = fortran_info(LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, *jobu, *jobvt, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt, work, *lwork));
}
void dpocon(const char* uplo, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info) {
    *info = fortran_info(LAPACKE_dpocon_work(LAPACK_COL_MAJOR, *uplo, *n, a, *lda, *anorm, rcond, work, iwork));
}
void dpotrf(const char* uplo, const int* n, double* a, const int* lda, int* info) {
    *info = fortran_info(LAPACKE_dpotrf_work(LAPACK_COL_MAJOR, *uplo, *n, a, *lda));
}
void dpotrs(const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda, double* b, const int* ldb, int* info) {
    *info = fortran_info(LAPACKE_dpotrs_work(LAPACK_COL_MAJOR, *uplo, *n, *nrhs, a, *lda, b, *ldb));
}
void dsyev(const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* w, double* work, const int* lwork, int* info) {
    *info = fortran_info(LAPACKE_dsyev_work(LAPACK_COL_MAJOR, *jobz, *uplo, *n, a, *lda, w, work, *lwork));
}
void dlarfg(const int* n, double* alpha, double* x, const int* incx, double* tau) {
    if (LAPACKE_dlarfg_work(*n, alpha, x, *incx, tau) != 0) throw std::runtime_error("LAPACK diffuse reflector failed");
}
void dgemm(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, const double* beta, double* c, const int* ldc) {
    cblas_dgemm(CblasColMajor, transpose(*transa), transpose(*transb), *m, *n, *k, *alpha, a, *lda, b, *ldb, *beta, c, *ldc);
}
const char* provider_name() noexcept { return "OpenBLAS"; }
const char* provider_version() noexcept { return openblas_get_config(); }
}}
