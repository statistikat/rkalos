#ifndef KALOS_DETAIL_LAPACK_HPP
#define KALOS_DETAIL_LAPACK_HPP
// Internal LP64 provider boundary. Exactly one implementation is linked.
// Owned column-major storage, caller-managed workspaces, Fortran INFO numbering.
namespace kalos { namespace lapack {
void dgelsy(const int* m, const int* n, const int* nrhs, double* a, const int* lda, double* b, const int* ldb, int* jpvt, const double* rcond, int* rank, double* work, const int* lwork, int* info);
void dgeev(const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda, double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr, double* work, const int* lwork, int* info);
void dgeqrf(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info);
void dgesvd(const char* jobu, const char* jobvt, const int* m, const int* n, double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt, double* work, const int* lwork, int* info);
void dpocon(const char* uplo, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info);
void dpotrf(const char* uplo, const int* n, double* a, const int* lda, int* info);
void dpotrs(const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda, double* b, const int* ldb, int* info);
void dsyev(const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* w, double* work, const int* lwork, int* info);
void dlarfg(const int* n, double* alpha, double* x, const int* incx, double* tau);
void dgemm(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, const double* beta, double* c, const int* ldc);
const char* provider_name() noexcept;
const char* provider_version() noexcept;
}}
#endif
