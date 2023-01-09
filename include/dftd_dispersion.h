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

/**
 * D4(EEQ)-ATM implementation
 */
#pragma once

#include "dftd_cutoff.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd4 {

class dparam {
 public:
  double s6;
  double s8;
  double s10;
  double s9;
  double a1;
  double a2;
  int alp;
};

/**
 * @brief Wrapper to handle the evaluation of dispersion energy and derivatives.
 * 
 * @param mol Molecular geometry.
 * @param charge Molecular charge.
 * @param par DFT-D4 parameters.
 * @param cutoff Real-space cutoffs for CN and dispersion.
 * @param energy Dispersion energy.
 * @param GRAD Dispersion gradient.
 * @return Exit status.
 */
extern int get_dispersion(
  const TMolecule &mol,
  const int charge,
  const dparam &par,
  TCutoff cutoff,
  double &energy,
  double *GRAD
);

// Generic wrappers for two- and three-body dispersion

extern int get_dispersion2(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient,
  bool lgrad = false
);

extern int get_dispersion3(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient,
  bool lgrad = false
);

// Two-body dispersion

extern int get_dispersion2_energy(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  TVector<double>& energy
);

extern int get_dispersion2_derivs(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient
);

// Three-body (ATM) dispersion

extern int get_atm_dispersion(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient,
  bool lgrad = false
); 

extern int get_atm_dispersion_energy(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double>& c6,
  TVector<double>& energy
);

extern int get_atm_dispersion_derivs(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient
);

extern double triple_scale(int ii, int jj, int kk);

extern inline double fdmpr_bj(const int n, const double r, const double c);
extern inline double fdmprdr_bj(const int n, const double r, const double c);

}  // namespace dftd4
