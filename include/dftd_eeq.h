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
 * Electronegativity equilibration (EEQ) model for DFT-D4
 */

#pragma once

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd4 {

extern int get_charges(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const int charge,
  const double cutoff,
  TVector<double> &q,
  TMatrix<double> &dqdr,
  bool lgrad
);

extern int get_vrhs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int &charge,
  const TVector<double> &cn,
  TVector<double> &Xvec,
  TVector<double> &dXvec,
  bool lgrad
);

extern int get_amat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TMatrix<double> &Amat
);

extern int get_damat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const TVector<double> &q,
  const TMatrix<double> &Amat,
  TMatrix<double> &dAmat,
  TMatrix<double> &atrace
);

extern int eeq_chrgeq(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const int &charge,
  const TVector<double> &cn,
  TVector<double> &qvec,
  TMatrix<double> &dcndr,
  TMatrix<double> &dqdr,
  bool lgrad = false,
  bool lverbose = false
);

} // namespace dftd4
