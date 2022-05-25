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

/**
 * D4(EEQ)-ATM implementation
 */
#pragma once

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd {

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

extern int d4dim(TMolecule& mol);

extern int d4(TMolecule& mol, int ndim, double wf, double g_a, double g_c,
              TVector<double>& cn, TVector<double>& gw, TMatrix<double>& c6ref);

extern int edisp(TMolecule& mol, TMatrix<double>& dist, int ndim,
                 TVector<double>& q, dparam& par, double g_a, double g_c,
                 TVector<double>& gw, TMatrix<double>& c6ref, bool lmbd,
                 double& energy);

extern int dispgrad(TMolecule& mol, TMatrix<double>& dist, int ndim,
                    TVector<double>& q, TMatrix<double>& dqdr,
                    TVector<double>& cn, TMatrix<double>& dcndr, dparam& par,
                    double wf, double g_a, double g_c, TMatrix<double>& c6ref,
                    bool lmbd, double& energy, TMatrix<double>& gradient);

extern int apprabc(TMolecule& mol, TMatrix<double>& dist, int ndim,
                   dparam& par, TVector<double>& c6ab, double& energy);

extern int dabcappr(TMolecule& mol, TMatrix<double>& dist, int ndim,
                    dparam& par, TVector<double>& gw, TVector<double>& dgw,
                    TMatrix<double>& c6ref, TVector<double>& dc6dr,
                    TVector<double>& dc6dcn, double& energy);

}  // namespace dftd
