/* This file is part of cpp-d4.
 *
 * Copyright (C) 2019 Sebastian Ehlert, Marvin Friede
 *
 * cpp-d4 is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cpp-d4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with cpp-d4.  If not, see <https://www.gnu.org/licenses/>.
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
  if (A.rows == C.N && A.cols == V.N) {
    if (Transpose) {
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
    } else {
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
};

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
      C.rows == 0)
    exit(EXIT_FAILURE);

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
};

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
  if (info != 0) { return EXIT_FAILURE; }

  // Inverse of an LU-factored general matrix
  info = LAPACKE_dgetri(
    LAPACK_ROW_MAJOR, (lapack_int)a.rows, a.p, (lapack_int)a.cols, ipiv
  );
  if (info != 0) { return EXIT_FAILURE; }

  delete[] ipiv;

  return EXIT_SUCCESS;
};

} // namespace dftd4
