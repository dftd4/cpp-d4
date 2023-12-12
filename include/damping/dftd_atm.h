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
 * @brief Three-body (ATM) dispersion
 */

#pragma once

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd4 {

extern int get_atm_dispersion(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double> &c6,
  const TMatrix<double> &dc6dcn,
  const TMatrix<double> &dc6dq,
  TVector<double> &energy,
  TVector<double> &dEdcn,
  TVector<double> &dEdq,
  TVector<double> &gradient,
  bool lgrad = false
);

extern int get_atm_dispersion_energy(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double> &c6,
  TVector<double> &energy
);

extern int get_atm_dispersion_derivs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double> &c6,
  const TMatrix<double> &dc6dcn,
  const TMatrix<double> &dc6dq,
  TVector<double> &energy,
  TVector<double> &dEdcn,
  TVector<double> &dEdq,
  TVector<double> &gradient
);

extern double triple_scale(int ii, int jj, int kk);

} // namespace dftd4
