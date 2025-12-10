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

#include <cstdlib>
#include <cstring>
#include <iostream>

namespace dftd4 {

// Define a vector
template <class T> class TVector {
  public:
    int N;           // Dimension of the vector
    int ElementSize; // Size of each element in the vector
    T *p;            // The pointer to the vector

    TVector() {
      N = 0;
      p = nullptr;
      ElementSize = sizeof(T);
    }
    ~TVector() {
      if (p != nullptr) Delete();
    }
    void NewVector(int VectorLength) {
      if (VectorLength < 0) { std::exit(EXIT_FAILURE); }
      if (p != nullptr && N == VectorLength) {
        Init();
      } else {
        Delete();
        if (VectorLength == 0) { return; }
        // get new memory
        p = new T[VectorLength];
        if (!p) { std::exit(EXIT_FAILURE); }
        N = VectorLength;
        Init();
      }
    }

    // alias for NewVector
    void New(int VectorLength) { return NewVector(VectorLength); }
    void NewVec(int VectorLength) { return NewVector(VectorLength); }

    void Delete() {
      if (p != nullptr && N != 0) { delete[] p; }
      p = nullptr;
      N = 0;
    }
    void DelVec(void) { return Delete(); }

    void CopyVec(const TVector &v) {
      long int mem;
      if (N != v.N) {
        Delete();
        New(v.N);
      }
      if (v.N == 0) return;
      mem = (long int)N * ElementSize;
      std::memcpy(p, v.p, mem);
    }

    void Init() {
      if (p != nullptr) {
        long int mem = (long int)N * ElementSize;
        std::memset(p, 0, mem);
      }
    }

    void Print(const char name[]) {
      printf("Vector printed: %s (%d)\n", name, N);
      for (int i = 0; i < N; i++) {
        printf("%+23.15e\n", p[i]);
      }
      printf("\n");
    }
    void PrintInt(const char name[]) {
      printf("Vector printed: %s (%d)\n", name, N);
      for (int i = 0; i < N; i++) {
        printf("%d\n", p[i]);
      }
      printf("\n");
    }

    inline T &operator()(int i) { return p[i]; }
    inline const T &operator()(int i) const { return p[i]; }
    inline T &operator[](int i) { return p[i]; }
    inline const T &operator[](int i) const { return p[i]; }

    // Max element
    T Max() const {
      T result = p[0];
      for (int i = 0; i < N; ++i)
        if (result < p[i]) result = p[i];
      return result;
    };
};

// Define a normal matrix
template <class T> class TMatrix {
  public:
    int rows, cols;  // dimensions
    int ElementSize; // Size of elements in matrix
    T *p;            // pointer to dynamic memory

    TMatrix() {
      cols = 0;
      rows = 0;
      p = nullptr;
      ElementSize = sizeof(T);
    }
    ~TMatrix() {
      if (p != nullptr) Delete();
    }

    void NewMatrix(int r, int c) {
      if (r < 0 || c < 0) std::exit(EXIT_FAILURE);
      if (p != nullptr && r == rows && c == cols) {
        Init();
      } else {
        long int mem = (long int)r * (long int)c;
        if (p != nullptr) Delete(); // Eventually delete old matrix

        if (mem == 0) return; // don't touch pointer if no memory is allocated

        p = new T[mem];
        if (!p) std::exit(EXIT_FAILURE);
        rows = r;
        cols = c;
        Init();
      }
      return;
    }

    void NewMatrix(const TMatrix &v) { NewMatrix(v.rows, v.cols); }

    // alias for NewMatrix
    void New(int r, int c) { return NewMatrix(r, c); }
    void NewMat(int r, int c) { return NewMatrix(r, c); }

    void Delete() {
      if (p != nullptr && rows * cols != 0) { delete[] p; }
      rows = 0;
      cols = 0;
      p = nullptr;
    }
    void DelMat(void) { return Delete(); }

    void Init() {
      long int mem;
      if (p != nullptr) {
        mem = (long int)cols * (long int)rows * ElementSize;
        std::memset(p, 0, mem);
      }
    }

    void Transpose() {
      T x;
      int i, j;

      if (p != nullptr) {
        if (rows == cols) {
          for (i = 0; i < rows; i++) {
            for (j = 0; j < i; j++) {
              x = p[i * cols + j];
              p[i * cols + j] = p[j * cols + i];
              p[j * cols + i] = x;
            } // j
          } // i
        } // if NxN
        else {
          // for non-square matrix, we need an additional copy
          TMatrix<T> temp;
          temp.CopyMat(*this);
          NewMatrix(cols, rows);
          for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
              p[i * cols + j] = temp.p[j * cols + i];
            } // j
          } // i
        }
      } // if data is loaded
    } // for NxN matrices transpose elements

    void CopyMat(const TMatrix &m) {
      long int mem;

      if ((m.rows != rows) || (m.cols != cols)) {
        Delete();
        New(m.rows, m.cols);
      }
      mem = (long int)rows * (long int)cols * ElementSize;
      if (mem == 0) return;
      std::memcpy(p, m.p, mem);
    }

    void Print(const char name[] = "unknown") {
      printf("Matrix printed: %s (%d, %d)\n", name, rows, cols);
      for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
          printf("%+23.15e", p[i * cols + j]);
          if (j == cols - 1) {
            printf("\n");
          } else {
            printf(" ");
          }
        }
      }
      printf("\n");
    }

    inline T &operator()(int i, int j) { return p[i * cols + j]; }
    inline const T &operator()(int i, int j) const { return p[i * cols + j]; }
    inline T *operator[](int i) { return p + i * cols; }
};

typedef TVector<int> TIVector;
typedef TVector<double> TRVector;

} // namespace dftd4
