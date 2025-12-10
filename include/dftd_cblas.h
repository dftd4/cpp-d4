/*
 * This file contains code adapted from the ORCA quantum chemistry program.
 * ORCA is developed by the group of Prof. Frank Neese at the
 * Max-Planck-Institut für Kohlenforschung, Mülheim an der Ruhr and FAccTs GmbH.
 * ORCA is licensed by the Max-Planck-Institut für Kohlenforschung and FAccTs
 * GmbH.
 *
 * The inclusion of ORCA code in this file has been done with the explicit
 * permission of the ORCA developers.
 *
 * For reuse or licensing of this code, please contact the ORCA team at the
 * Max-Planck-Institut für Kohlenforschung (https://orcaforum.kofo.mpg.de/) or
 * FAccTs GmbH (https://www.faccts.de/).
 */
#pragma once

#include "dftd_matrix.h"

#include "cblas.h"
#include "lapacke.h"

namespace dftd4 {

/**
 * @brief General matrix vector multiplication (`C = alpha * A * V + C`).
 *
 * @param C Result vector C. Modified in-place.
 * @param A Matrix A.
 * @param V Vector V.
 * @param Transpose Specifies whether to transpose matrix A.
 * @param alpha Scaling factor for the product of matrix A and vector X.
 * @return Exit code
 */
inline int BLAS_Add_Mat_x_Vec(
  TVector<double> &C,
  TMatrix<double> &A,
  TVector<double> &V,
  bool Transpose,
  double alpha
) {
  if (Transpose) {
    if (A.cols == C.N && A.rows == V.N) {
      cblas_dgemv(
        CblasRowMajor,
        CblasTrans,
        A.rows,
        A.cols,
        alpha,
        A.p,
        A.cols,
        V.p,
        1,
        1.0,
        C.p,
        1
      );
      return EXIT_SUCCESS;
    };
  } else {
    if (A.rows == C.N && A.cols == V.N) {
      cblas_dgemv(
        CblasRowMajor,
        CblasNoTrans,
        A.rows,
        A.cols,
        alpha,
        A.p,
        A.cols,
        V.p,
        1,
        1.0,
        C.p,
        1
      );
      return EXIT_SUCCESS;
    };
  };

  return EXIT_FAILURE;
}

/**
 * @brief General matrix-matrix multiplication (`C = alpha * A * B + C`).
 *
 * @param C Result matrix C. Modified in-place.
 * @param A Matrix A.
 * @param B Matrix B.
 * @param TransposeA Specifies whether to transpose matrix A.
 * @param TransposeB Specifies whether to transpose matrix B.
 * @param alpha Scaling factor for the product of matrix A and matrix B.
 * @return Exit code.
 */
inline int BLAS_Add_Mat_x_Mat(
  TMatrix<double> &C,
  const TMatrix<double> &A,
  const TMatrix<double> &B,
  const bool TransposeA,
  const bool TransposeB,
  const double alpha
) {
  // check for size 0 matrices
  if (A.cols == 0 || A.rows == 0 || B.cols == 0 || B.rows == 0 || C.cols == 0 ||
      C.rows == 0) {
    exit(EXIT_FAILURE);
  };

  // check for transpositions
  if (!TransposeA) {
    if (!TransposeB) {
      // check dimensions
      if (A.cols != B.rows || A.rows != C.rows || B.cols != C.cols) {
        exit(EXIT_FAILURE);
      };
      cblas_dgemm(
        CblasRowMajor,
        CblasNoTrans,
        CblasNoTrans,
        C.rows,
        C.cols,
        A.cols,
        alpha,
        A.p,
        A.cols,
        B.p,
        B.cols,
        1.0,
        C.p,
        C.cols
      );
    } // B not transposed
    else {
      // check dimensions for C=A*BT
      if (A.cols != B.cols || A.rows != C.rows || B.rows != C.cols) {
        exit(EXIT_FAILURE);
      };
      // B is transposed, A not
      cblas_dgemm(
        CblasRowMajor,
        CblasNoTrans,
        CblasTrans,
        C.rows,
        C.cols,
        A.cols,
        alpha,
        A.p,
        A.cols,
        B.p,
        B.cols,
        1.0,
        C.p,
        C.cols
      );
    }; // B transposed
  } // A not transposed
  else {
    if (!TransposeB) {
      // check dimensions for C=AT*B
      if (A.rows != B.rows || A.cols != C.rows || B.cols != C.cols) {
        exit(EXIT_FAILURE);
      };
      // A is transposed and B not
      cblas_dgemm(
        CblasRowMajor,
        CblasTrans,
        CblasNoTrans,
        C.rows,
        C.cols,
        A.rows,
        alpha,
        A.p,
        A.cols,
        B.p,
        B.cols,
        1.0,
        C.p,
        C.cols
      );
    } // B not transposed
    else {
      // check dimensions for C=AT*BT
      if (A.rows != B.cols || A.cols != C.rows || B.rows != C.cols) {
        exit(EXIT_FAILURE);
      };
      // both are transposed
      cblas_dgemm(
        CblasRowMajor,
        CblasTrans,
        CblasTrans,
        C.rows,
        C.cols,
        A.rows,
        alpha,
        A.p,
        A.cols,
        B.p,
        B.cols,
        1.0,
        C.p,
        C.cols
      );
    }; // B transposed
  };
  return EXIT_SUCCESS;
}

/**
 * @brief Compute inverse of a matrix using LU decomposition.
 *
 * @param a Matrix a.
 * @return Exit code.
 */
inline int BLAS_InvertMatrix(TMatrix<double> &a) {
  if (a.rows != a.cols) { return EXIT_FAILURE; }

  lapack_int info;
  lapack_int *ipiv = new lapack_int[a.rows];

  // LU factorization of a general m-by-n matrix
  info = LAPACKE_dgetrf(
    LAPACK_ROW_MAJOR,
    (lapack_int)a.rows,
    (lapack_int)a.cols,
    a.p,
    (lapack_int)a.cols,
    ipiv
  );
  if (info != 0) {
    delete[] ipiv;
    return EXIT_FAILURE;
  }

  // Inverse of an LU-factored general matrix
  info = LAPACKE_dgetri(
    LAPACK_ROW_MAJOR, (lapack_int)a.rows, a.p, (lapack_int)a.cols, ipiv
  );
  delete[] ipiv;

  if (info != 0) { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}

/**
 * @brief Solve a symmetric linear system A * X = B for X.
 *
 * This routine factorizes a symmetric matrix A using Bunch-Kaufman
 * factorization and solves for the right-hand side vector B. The matrix A is
 * overwritten by its factorization. The solution overwrites B.
 *
 * @param A Symmetric matrix of size (m x m). Overwritten by the factorization.
 * @param B Right-hand side vector of size m. Overwritten by the solution.
 * @return int Returns EXIT_SUCCESS (0) on success, EXIT_FAILURE (1) on error.
 */
inline int BLAS_SolveSymmetric(
  TMatrix<double> &A, // symmetric matrix
  TVector<double> &B  // RHS vector (becomes solution)
) {
  const lapack_int m = A.rows;
  const lapack_int nrhs = 1;

  if (A.cols != m || B.N != m) {
    fprintf(stderr, "BLAS_SolveSymmetric error: dimension mismatch\n");
    return EXIT_FAILURE;
  }

  lapack_int info;
  lapack_int *ipiv = new lapack_int[A.rows];

  // Factorization
  info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'L', m, A.p, m, ipiv);
  if (info != 0) {
    delete[] ipiv;
    fprintf(stderr, "dsytrf failed: info=%d\n", (int)info);
    return EXIT_FAILURE;
  }

  // Solve for all RHS columns
  info =
    LAPACKE_dsytrs(LAPACK_ROW_MAJOR, 'L', m, nrhs, A.p, m, ipiv, B.p, nrhs);
  delete[] ipiv;

  if (info != 0) {
    fprintf(stderr, "dsytrs failed: info=%d\n", (int)info);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

} // namespace dftd4
