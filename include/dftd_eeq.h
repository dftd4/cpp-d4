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

/**
 * Electronegativity equilibration (EEQ) model for DFT-D4
 */

#pragma once

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd {

extern int eeq_chrgeq(TMolecule& mol, int& charge, TMatrix<double>& dist,
                      TVector<double>& cn, TMatrix<double>& dcndr,
                      TVector<double>& q, TMatrix<double>& dqdr, double& energy,
                      TMatrix<double>& gradient, bool lverbose = false,
                      bool lgrad = false, bool lcpq = false);

}
