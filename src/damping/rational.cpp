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

#include <cmath>

#include "dftd_cblas.h"
#include "dftd_eeq.h"
#include "dftd_dispersion.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_ncoord.h"
#include "dftd_parameters.h"
#include "damping/atm.h"
#include "damping/rational.h"


namespace dftd4 {

inline double fdmpr_bj(const int n, const double r, const double c) {
  return 1.0 / (pow(r, n) + pow(c, n));
}
inline double fdmprdr_bj(const int n, const double r, const double c) {
  return -n * pow(r, n - 1) * pow(fdmpr_bj(n, r, c), 2);
}

int get_dispersion2(
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
  bool lgrad /*= false*/
) {
  int info{0};

  if (lgrad) {
    info = get_dispersion2_derivs(
      mol, dist, cutoff, par, c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient
    );
  } else {
    info = get_dispersion2_energy(mol, dist, cutoff, par, c6, energy);
  }

  return info;
}

int get_dispersion2_energy(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  TVector<double>& energy
) {
  int izp{0}, jzp{0};
  double r{0.0}, r0ij{0.0}, r4r2ij{0.0}, c6ij{0.0};
  double t6{0.0}, t8{0.0}, t10{0.0};
  double e{0.0}, edisp{0.0};

  for (int iat = 0; iat != mol.NAtoms; iat++) {
    izp = mol.at(iat);
    for (int jat = 0; jat != iat; jat++) {
      jzp = mol.at(jat);
      r = dist(iat, jat);
      if (r > cutoff) continue;
      
      r4r2ij = 3.0 * r4r2[izp] * r4r2[jzp];
      r0ij = par.a1 * sqrt(r4r2ij) + par.a2;
      c6ij = c6(iat, jat);

      t6 = fdmpr_bj(6, r, r0ij);
      t8 = fdmpr_bj(8, r, r0ij);

      edisp = par.s6 * t6 + par.s8 * r4r2ij * t8;

      if (par.s10 != 0.0) {
        t10 = fdmpr_bj(10, r, r0ij);
        edisp += par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * t10;
      }

      e = -c6ij*edisp * 0.5;
      energy(iat) += e;
      energy(jat) += e;
    }
  }

  return EXIT_SUCCESS;
}

int get_dispersion2_derivs(
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
) {
  int izp{0}, jzp{0};
  double r{0.0}, r0ij{0.0}, r4r2ij{0.0}, c6ij{0.0};
  double t6{0.0}, t8{0.0}, t10{0.0};
  double dt6{0.0}, dt8{0.0}, dt10{0.0};
  double e{0.0}, edisp{0.0}, gdisp{0.0};
  double x{0.0}, y{0.0}, z{0.0};
  double dgx{0.0}, dgy{0.0}, dgz{0.0};

  for (int iat = 0; iat != mol.NAtoms; iat++) {
    izp = mol.at(iat);
    for (int jat = 0; jat != iat; jat++) {
      jzp = mol.at(jat);
      r = dist(iat, jat);
      if (r > cutoff) continue;
      
      r4r2ij = 3.0 * r4r2[izp] * r4r2[jzp];
      r0ij = par.a1 * sqrt(r4r2ij) + par.a2;
      c6ij = c6(iat, jat);

      t6 = fdmpr_bj(6, r, r0ij);
      t8 = fdmpr_bj(8, r, r0ij);

      dt6 = -6 * pow(r, 5) * pow(t6, 2);
      dt8 = -8 * pow(r, 7) * pow(t8, 2);

      edisp = par.s6 * t6 + par.s8 * r4r2ij * t8;
      gdisp = par.s6 * dt6 + par.s8 * r4r2ij * dt8;

      if (par.s10 != 0.0) {
        t10 = fdmpr_bj(10, r, r0ij);
        dt10 = -10 * pow(r, 9) * pow(t10, 2);
        edisp += par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * t10;
        gdisp += par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * dt10;
      }

      e = -c6ij*edisp * 0.5;
      energy(iat) += e;
      energy(jat) += e;

      x = (mol.xyz(iat, 0) - mol.xyz(jat, 0)) / r;
      y = (mol.xyz(iat, 1) - mol.xyz(jat, 1)) / r;
      z = (mol.xyz(iat, 2) - mol.xyz(jat, 2)) / r;
      dgx = -c6ij * gdisp * x;
      dgy = -c6ij * gdisp * y;
      dgz = -c6ij * gdisp * z;

      dEdcn(iat) -= dc6dcn(iat, jat) * edisp;
      dEdcn(jat) -= dc6dcn(jat, iat) * edisp;

      dEdq(iat) -= dc6dq(iat, jat) * edisp;
      dEdq(jat) -= dc6dq(jat, iat) * edisp;

      gradient(3*iat) += dgx;
      gradient(3*iat + 1) += dgy;
      gradient(3*iat + 2) += dgz;
      gradient(3*jat) -= dgx;
      gradient(3*jat + 1) -= dgy;
      gradient(3*jat + 2) -= dgz;
    }
  }
  return EXIT_SUCCESS;
}


int get_dispersion3(
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
  bool lgrad /*= false*/
) {
  return get_atm_dispersion(
    mol, dist, cutoff, par.s9, par.a1, par.a2, par.alp, c6, dc6dcn, dc6dq,
    energy, dEdcn, dEdq, gradient, lgrad
  );
}

}  // namespace dftd4
