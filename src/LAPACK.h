#ifndef qpOASES_LAPACK_H
#define qpOASES_LAPACK_H

#ifdef __USE_SINGLE_PRECISION__

	/** Macro for calling level 3 BLAS operation in single precision. */
	#define GEMM sgemm_
	/** Macro for calling level 3 BLAS operation in single precision. */
	#define SYR ssyr_
	/** Macro for calling level 3 BLAS operation in single precision. */
	#define SYR2 ssyr2_
	/** Macro for calling level 3 BLAS operation in single precision. */
	#define POTRF spotrf_
  /** Macro for calling level 3 BLAS operation in single precision. */
  //#define GEQRF sgeqrf_
  /** Macro for calling level 3 BLAS operation in single precision. */
  //#define ORMQR sormqr_
  /** Macro for calling level 3 BLAS operation in single precision. */
  #define TRTRS strtrs_
  /** Macro for calling level 3 BLAS operation in single precision. */
  #define TRCON strcon_

#else

	/** Macro for calling level 3 BLAS operation in double precision. */
	#define GEMM dgemm_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define SYR  dsyr_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define SYR2 dsyr2_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define POTRF dpotrf_
  /** Macro for calling level 3 BLAS operation in double precision. */
	//#define GEQRF dgeqrf_
	/** Macro for calling level 3 BLAS operation in double precision. */
	//#define ORMQR dormqr_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define TRTRS dtrtrs_
	/** Macro for calling level 3 BLAS operation in double precision. */
	#define TRCON dtrcon_

#endif /* __USE_SINGLE_PRECISION__ */


extern "C"
{
	/** Performs one of the matrix-matrix operation in double precision. */
	void dgemm_ ( const char*, const char*, const unsigned long*, const unsigned long*, const unsigned long*,
			const double*, const double*, const unsigned long*, const double*, const unsigned long*,
			const double*, double*, const unsigned long* );
	/** Performs one of the matrix-matrix operation in single precision. */
	void sgemm_ ( const char*, const char*, const unsigned long*, const unsigned long*, const unsigned long*,
			const float*, const float*, const unsigned long*, const float*, const unsigned long*,
			const float*, float*, const unsigned long* );

	/** Performs a symmetric rank 1 operation in double precision. */
	void dsyr_ ( const char *, const unsigned long *, const double *, const double *,
				 const unsigned long *, double *, const unsigned long *);
	/** Performs a symmetric rank 1 operation in single precision. */
	void ssyr_ ( const char *, const unsigned long *, const float *, const float *,
				 const unsigned long *, float *, const unsigned long *);

	/** Performs a symmetric rank 2 operation in double precision. */
	void dsyr2_ ( const char *, const unsigned long *, const double *, const double *,
				  const unsigned long *, const double *, const unsigned long *, double *, const unsigned long *);
	/** Performs a symmetric rank 2 operation in single precision. */
	void ssyr2_ ( const char *, const unsigned long *, const float *, const float *,
				  const unsigned long *, const float *, const unsigned long *, float *, const unsigned long *);

	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in double precision. */
	void dpotrf_ ( const char *, const unsigned long *, double *, const unsigned long *, long * );
	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in single precision. */
	void spotrf_ ( const char *, const unsigned long *, float *, const unsigned long *, long * );

	/** Compute a QR factorization of a real M-by-N matrix A in double precision */
	//void dgeqrf_(	const unsigned long *M, const unsigned long *N, double *A, const unsigned long *LDA,
					//double *TAU, double *WORK, const unsigned long *LWORK, int *INFO );
	/** Compute a QR factorization of a real M-by-N matrix A in single precision */
	//void sgeqrf_(	const unsigned long *M, const unsigned long *N, float *A, const unsigned long *LDA,
					//float *TAU, float *WORK, const unsigned long *LWORK, int *INFO );

	/** Multiply C with orthogonal matrix Q**T as returned by geqrf (double precision) */
	//void dormqr_(	const char *SIDE, const char *TRANS, const unsigned long *M, const unsigned long *N, const unsigned long *K,
					//double *A, const unsigned long *LDA, double *TAU, double *C, const unsigned long *LDC,
					//double *WORK, const unsigned long *LWORK, int *INFO );
	/** Multiply C with orthogonal matrix Q**T as returned by geqrf (single precision) */
	//void sormqr_(	const char *SIDE, const char *TRANS, const unsigned long *M, const unsigned long *N, const unsigned long *K,
					//float *A, const unsigned long *LDA, float *TAU, float *C, const unsigned long *LDC,
					//float *WORK, const unsigned long *LWORK, int *INFO );

	/** Solve a triangular system (double precision) */
	void dtrtrs_(	const char *UPLO, const char *TRANS, const char *DIAG, const unsigned long *N, const unsigned long *NRHS,
					double *A, const unsigned long *LDA, double *B, const unsigned long *LDB, long *INFO );
	/** Solve a triangular system (single precision) */
	void strtrs_(	const char *UPLO, const char *TRANS, const char *DIAG, const unsigned long *N, const unsigned long *NRHS,
					float *A, const unsigned long *LDA, float *B, const unsigned long *LDB, long *INFO );

	/** Estimate the reciprocal of the condition number of a triangular matrix in double precision */
	void dtrcon_(	const char *NORM, const char *UPLO, const char *DIAG, const unsigned long *N, double *A, const unsigned long *LDA,
					double *RCOND, double *WORK, const unsigned long *IWORK, long *INFO );
	/** Estimate the reciprocal of the condition number of a triangular matrix in single precision */
	void strcon_(	const char *NORM, const char *UPLO, const char *DIAG, const unsigned long *N, float *A, const unsigned long *LDA,
					float *RCOND, float *WORK, const unsigned long *IWORK, long *INFO );
}

#endif // qpOASES_LAPACK_H
