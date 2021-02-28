/* This file is part of cpp-d4.
 *
 * Copyright (C) 2019-2021 Sebastian Ehlert
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

/*
 * Definition of the coordination numbers used in DFT-D4 and the
 * electronegativity equilibration (EEQ) model
 */

#pragma once

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd {
// Calculate all distance pairs and store in matrix
extern int calc_distances(TMolecule& mol, TMatrix<double>& dist);

// Calculate error function coordination number
extern int ncoord_erf(TMolecule& mol, TMatrix<double>& dist,
                      TVector<double>& cn, double thr = 1600.0);

// Calculate covalent coordination number
extern int ncoord_d4(TMolecule& mol, TMatrix<double>& dist,
                     TVector<double>& cn, double thr = 1600.0);

// Calculate error function coordination number and derivative
// w.r.t. nuclear coordinates
extern int dncoord_erf(TMolecule& mol, TMatrix<double>& dist,
                       TVector<double>& cn, TMatrix<double>& dcndr,
                       double thr = 1600.0);

// Calculate covalent coordination number and derivative
// w.r.t. nuclear coordinates
extern int dncoord_d4(TMolecule& mol, TMatrix<double>& dist,
                      TVector<double>& cn, TMatrix<double>& dcndr,
                      double thr = 1600.0);

};  // namespace dftd
