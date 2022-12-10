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

extern
int d4dim(const TMolecule& mol);

extern
int d4(const TMolecule& mol, int ndim, double wf, double g_a, double g_c,
       const TVector<double>& cn, TVector<double>& gw, TMatrix<double>& c6ref);

extern
int edisp(const TMolecule& mol, const TMatrix<double>& dist, const dparam& par, 
          int ndim, TVector<double>& q,double g_a, double g_c,
          TVector<double>& gw, TMatrix<double>& c6ref, bool lmbd,
          double& energy);

extern
int dispgrad(const TMolecule& mol, const TMatrix<double>& dist,
             const dparam& par, int ndim, const TVector<double>& q, 
             TMatrix<double>& dqdr, TVector<double>& cn, 
             TMatrix<double>& dcndr, double wf, double g_a, double g_c, 
             TMatrix<double>& c6ref, bool lmbd, double& energy, 
             TMatrix<double>& gradient);

extern
int apprabc(const TMolecule& mol, const TMatrix<double>& dist,
            const dparam& par, int ndim, TVector<double>& c6ab,
            double& energy);

extern
int dabcappr(const TMolecule& mol, const TMatrix<double>& dist, 
             const dparam& par, int ndim, TVector<double>& gw, 
             TVector<double>& dgw, TMatrix<double>& c6ref, 
             TVector<double>& dc6dr, TVector<double>& dc6dcn, double& energy);

extern
int DFTVDW_D4(const TMolecule &mol, const dparam &par, const int &charge,
              double &energy, double *GRAD);

}  // namespace dftd
