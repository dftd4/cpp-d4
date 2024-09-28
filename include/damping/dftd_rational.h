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

#include "dftd_dispersion.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd4 {

extern inline double fdmpr_bj(int n, double r, double c);
extern inline double fdmprdr_bj(int n, double r, double c);

// Generic wrappers for two- and three-body dispersion

extern int get_dispersion2(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  const dparam &par,
  const TMatrix<double> &c6,
  const TMatrix<double> &dc6dcn,
  const TMatrix<double> &dc6dq,
  TVector<double> &energy,
  TVector<double> &dEdcn,
  TVector<double> &dEdq,
  TVector<double> &gradient,
  bool lgrad = false
);

extern int get_dispersion3(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  const dparam &par,
  const TMatrix<double> &c6,
  const TMatrix<double> &dc6dcn,
  const TMatrix<double> &dc6dq,
  TVector<double> &energy,
  TVector<double> &dEdcn,
  TVector<double> &dEdq,
  TVector<double> &gradient,
  bool lgrad = false
);

// Two-body dispersion

extern int get_dispersion2_energy(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  const dparam &par,
  const TMatrix<double> &c6,
  TVector<double> &energy
);

extern int get_dispersion2_derivs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  const dparam &par,
  const TMatrix<double> &c6,
  const TMatrix<double> &dc6dcn,
  const TMatrix<double> &dc6dq,
  TVector<double> &energy,
  TVector<double> &dEdcn,
  TVector<double> &dEdq,
  TVector<double> &gradient
);

} // namespace dftd4
