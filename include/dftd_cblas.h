/* This file is part of cpp-d4.
 *
 * Copyright (C) 2019-2020 Sebastian Ehlert
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

extern "C" {
#include "cblas.h"
#include "lapacke.h"
}

namespace dftd {

inline int BLAS_Add_Mat_x_Vec(TVector<double>& C, TMatrix<double>& A,
                              TVector<double>& V, bool Transpose,
                              double alpha) {
  if (A.rows == A.cols) {
    if (Transpose) {
      cblas_dgemv(CblasRowMajor, CblasTrans, A.rows, A.cols, alpha, A.p,
                  A.cols, V.p, 1, 1.0, C.p, 1);
      return 0;
    } else {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, A.rows, A.cols, alpha, A.p,
                  A.cols, V.p, 1, 1.0, C.p, 1);
      return 0;
    };
  };

  return 0;
};

inline void BLAS_Add_Mat_x_Vec(double* C, const double* A, const double* B,
                               const int rows, const int cols,
                               const bool Transpose, const double alpha,
                               const double beta) {
  if (Transpose) {
    cblas_dgemv(CblasRowMajor, CblasTrans, rows, cols, alpha, A, cols, B, 1,
                beta, C, 1);
  } else {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rows, cols, alpha, A, cols, B, 1,
                beta, C, 1);
  };
};

inline double BLAS_Vec_x_Vec(const TVector<double>& A,
                             const TVector<double>& B) {
  if (A.N != B.N) exit(EXIT_FAILURE);
  return cblas_ddot(A.N, A.p, 1, B.p, 1);
};

/// C = alpha * A * B + C
inline int BLAS_Add_Mat_x_Mat(TMatrix<double>& C, const TMatrix<double>& A,
                              const TMatrix<double>& B, const bool TransposeA,
                              const bool TransposeB, const double alpha) {
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
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, C.rows, C.cols,
                  A.cols, alpha, A.p, A.cols, B.p, B.cols, 1.0, C.p, C.cols);
    }  // B not transposed
    else {
      // check dimensions for C=A*BT
      if (A.cols != B.cols || A.rows != C.rows || B.rows != C.cols) {
        exit(EXIT_FAILURE);
      };
      // B is transposed, A not
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, C.rows, C.cols,
                  A.cols, alpha, A.p, A.cols, B.p, B.cols, 1.0, C.p, C.cols);
    };  // B transposed
  }     // A not transposed
  else {
    if (!TransposeB) {
      // check dimensions for C=AT*B
      if (A.rows != B.rows || A.cols != C.rows || B.cols != C.cols) {
        exit(EXIT_FAILURE);
      };
      // A is transposed and B not
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C.rows, C.cols,
                  A.rows, alpha, A.p, A.cols, B.p, B.cols, 1.0, C.p, C.cols);
    }  // B not transposed
    else {
      // check dimensions for C=AT*BT
      if (A.rows != B.cols || A.cols != C.rows || B.rows != C.cols) {
        exit(EXIT_FAILURE);
      };
      // both are transposed
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, C.rows, C.cols, A.rows,
                  alpha, A.p, A.cols, B.p, B.cols, 1.0, C.p, C.cols);
    };  // B transposed
  };
  return EXIT_SUCCESS;
};

// Linear equation solver
inline int BLAS_LINEQ(TMatrix<double>& a, TMatrix<double>& b, int m) {
  lapack_int info, n, nrhs;
  lapack_int* ipiv;

  if (a.rows != a.cols) return EXIT_FAILURE;

  n = a.rows;
  nrhs = m;

  a.Transpose();
  b.Transpose();

  ipiv = new lapack_int[n];

  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, a.p, n, ipiv, b.p, n);

  delete[] ipiv;

  a.Transpose();
  b.Transpose();

  return info;
}

inline int BLAS_InvertMatrix(TMatrix<double>& a) {
  if (a.rows != a.cols) {
    return EXIT_FAILURE;
  }

  lapack_int info;
  lapack_int* ipiv = new lapack_int[a.rows];

  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, (lapack_int)a.rows,
                        (lapack_int)a.cols, a.p, (lapack_int)a.cols, ipiv);

  if (info != 0) {
    return EXIT_FAILURE;
  }

  info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, (lapack_int)a.rows, a.p,
                        (lapack_int)a.cols, ipiv);

  if (info != 0) {
    return EXIT_FAILURE;
  }

  delete[] ipiv;

  return EXIT_SUCCESS;
};

}  // namespace dftd
