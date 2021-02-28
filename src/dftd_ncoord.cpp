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
 * Definition of the coordination numbers used in DFT-D4 and the
 * electronegativity equilibration (EEQ) model
 *
 * This module works on a distance matrix to avoid recalculating
 * the distances every time.
 */

#include "dftd_ncoord.h"

#include <cmath>

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd {

// convert bohr (a.u.) to Ångström and back
static const double autoaa = 0.52917726;
static const double aatoau = 1.0 / autoaa;

// covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
// 188-197), values for metals decreased by 10 %
static const double rad[119]{
    0.0,  0.32, 0.46,                                            //  H,He
    1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,              //  Li-Ne
    1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,              //  Na-Ar
    1.76, 1.54,                                                  //  K,Ca
    1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,  //  Sc-Zn
    1.12, 1.09, 1.15, 1.10, 1.14, 1.17,                          //  Ga-Kr
    1.89, 1.67,                                                  //  Rb,Sr
    1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23,  //  Y-Cd
    1.28, 1.26, 1.26, 1.23, 1.32, 1.31,                          //  In-Xe
    2.09, 1.76,                                                  //  Cs,Ba
    1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51,                    //  La-Eu
    1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,                    //  Gd-Yb
    1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,  //  Lu-Hg
    1.30, 1.30, 1.36, 1.31, 1.38, 1.42,                          //  Tl-Rn
    2.01, 1.81,                                                  //  Fr,Ra
    1.67, 1.58, 1.52, 1.53, 1.54, 1.55, 1.49,                    //  Ac-Am
    1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58,                    //  Cm-No
    1.45, 1.41, 1.34, 1.29, 1.27, 1.21, 1.16, 1.15, 1.09, 1.22,  //  Lr-Cn
    1.22, 1.29, 1.46, 1.58, 1.48, 1.41                           //  Nh-Og
};

// pauling EN's
static const double en[119]{
    0.0, 2.20, 3.00,                                             //  H,He
    0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 4.50,              //  Li-Ne
    0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 3.50,              //  Na-Ar
    0.82, 1.00,                                                  //  K,Ca
    1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65,  //  Sc-Zn
    1.81, 2.01, 2.18, 2.55, 2.96, 3.00,                          //  Ga-Kr
    0.82, 0.95,                                                  //  Rb,Sr
    1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69,  //  Y-Cd
    1.78, 1.96, 2.05, 2.10, 2.66, 2.60,                          //  In-Xe
    0.79, 0.89,                                                  //  Cs,Ba
    1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18,                    //  La-Eu
    1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26,                    //  Gd-Yb
    1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00,  //  Lu-Hg
    1.62, 2.33, 2.02, 2.00, 2.20, 2.20,                          //  Tl-Rn
    // only dummies below
    1.50, 1.50,                                                  //  Fr,Ra
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                    //  Ac-Am
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                    //  Cm-No
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,  //  Rf-Cn
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50                           //  Nh-Og
};

static const double kn = 7.5;
static const double k4 = 4.10451;
static const double k5 = 19.08857;
static const double k6 = 2 * pow(11.28174, 2);
static const double hlfosqrtpi = 1.0 / 1.77245385091;

/*
 * Calculate all distance pairs and store in matrix
 */
int calc_distances(TMolecule& mol, TMatrix<double>& dist) {
  double rx = 0.0, ry = 0.0, rz = 0.0, tmp = 0.0;
  for (int i = 0; i != mol.NAtoms; i++) {
    dist(i, i) = 0.0;
    for (int j = 0; j != i; j++) {
      rx = mol.xyz(i, 0) - mol.xyz(j, 0);
      ry = mol.xyz(i, 1) - mol.xyz(j, 1);
      rz = mol.xyz(i, 2) - mol.xyz(j, 2);
      tmp = sqrt(rx * rx + ry * ry + rz * rz);
      dist(i, j) = tmp;
      dist(j, i) = tmp;
    }
  }

  return EXIT_SUCCESS;
}

/*
 * Calculate covalent coordination number
 */
int ncoord_d4(TMolecule& mol, TMatrix<double>& dist, TVector<double>& cn,
              double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double den = 0.0;
  double tmp = 0.0;

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != i; j++) {
      r = dist(i, j);
      rcovij = rad[mol.at(i)] + rad[mol.at(j)];
      rr = 3.0 * r / (4.0 * rcovij * aatoau);
      den =
          k4 *
          std::exp(-pow((abs(en[mol.at(i)] - en[mol.at(j)]) + k5), 2) / k6);
      tmp = den * 0.5 * (1.0 + std::erf(-kn * (rr - 1.0)));
      cn(i) += tmp;
      cn(j) += tmp;
    }
  }
  return EXIT_SUCCESS;
}

/*
 * Calculate covalent coordination number and derivative w.r.t. nuclear
 * coordinates
 */
int dncoord_d4(TMolecule& mol, TMatrix<double>& dist, TVector<double>& cn,
               TMatrix<double>& dcndr, double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0;
  double tmp = 0.0, dtmp = 0.0, den = 0.0;

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != i; j++) {
      r = dist(i, j);
      rx = (mol.xyz(j, 0) - mol.xyz(i, 0)) / r;
      ry = (mol.xyz(j, 1) - mol.xyz(i, 1)) / r;
      rz = (mol.xyz(j, 2) - mol.xyz(i, 2)) / r;
      rcovij = 4.0 / 3.0 * aatoau * (rad[mol.at(i)] + rad[mol.at(j)]);
      rr = r / rcovij;
      den =
          k4 *
          std::exp(-pow((abs(en[mol.at(i)] - en[mol.at(j)]) + k5), 2) / k6);
      tmp = den * 0.5 * (1.0 + std::erf(-kn * (rr - 1.0)));
      cn(i) += tmp;
      cn(j) += tmp;
      dtmp =
          -den * hlfosqrtpi * kn * std::exp(-pow(kn * (rr - 1.0), 2)) / rcovij;
      dcndr(j, 3 * j) += dtmp * rx;
      dcndr(j, 3 * j + 1) += dtmp * ry;
      dcndr(j, 3 * j + 2) += dtmp * rz;
      dcndr(j, 3 * i) = dtmp * rx;
      dcndr(j, 3 * i + 1) = dtmp * ry;
      dcndr(j, 3 * i + 2) = dtmp * rz;
      dcndr(i, 3 * j) = -dtmp * rx;
      dcndr(i, 3 * j + 1) = -dtmp * ry;
      dcndr(i, 3 * j + 2) = -dtmp * rz;
      dcndr(i, 3 * i) += -dtmp * rx;
      dcndr(i, 3 * i + 1) += -dtmp * ry;
      dcndr(i, 3 * i + 2) += -dtmp * rz;
    }
  }
  return EXIT_SUCCESS;
}

/*
 * Coordination number with error function based counting
 */
int ncoord_erf(TMolecule& mol, TMatrix<double>& dist, TVector<double>& cn,
               double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double tmp = 0.0;

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != i; j++) {
      r = dist(i, j);
      rcovij = rad[mol.at(i)] + rad[mol.at(j)];
      rr = 3.0 * r / (4.0 * rcovij * aatoau);
      tmp = 0.5 * (1.0 + std::erf(-kn * (rr - 1.0)));
      cn(i) += tmp;
      cn(j) += tmp;
    }
  }
  return EXIT_SUCCESS;
}

/*
 * Calculate error function coordination number and derivative w.r.t. nuclear
 * coordinates
 */
int dncoord_erf(TMolecule& mol, TMatrix<double>& dist, TVector<double>& cn,
                TMatrix<double>& dcndr, double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0;
  double tmp = 0.0, dtmp = 0.0;

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != i; j++) {
      r = dist(i, j);
      rx = (mol.xyz(j, 0) - mol.xyz(i, 0)) / r;
      ry = (mol.xyz(j, 1) - mol.xyz(i, 1)) / r;
      rz = (mol.xyz(j, 2) - mol.xyz(i, 2)) / r;
      rcovij = 4.0 / 3.0 * aatoau * (rad[mol.at(i)] + rad[mol.at(j)]);
      rr = r / rcovij;
      tmp = 0.5 * (1.0 + erf(-kn * (rr - 1.0)));
      cn(i) += tmp;
      cn(j) += tmp;
      dtmp = -hlfosqrtpi * kn * std::exp(-pow(kn * (rr - 1.0), 2)) / rcovij;
      dcndr(j, 3 * j) += dtmp * rx;
      dcndr(j, 3 * j + 1) += dtmp * ry;
      dcndr(j, 3 * j + 2) += dtmp * rz;
      dcndr(j, 3 * i) = dtmp * rx;
      dcndr(j, 3 * i + 1) = dtmp * ry;
      dcndr(j, 3 * i + 2) = dtmp * rz;
      dcndr(i, 3 * j) = -dtmp * rx;
      dcndr(i, 3 * j + 1) = -dtmp * ry;
      dcndr(i, 3 * j + 2) = -dtmp * rz;
      dcndr(i, 3 * i) += -dtmp * rx;
      dcndr(i, 3 * i + 1) += -dtmp * ry;
      dcndr(i, 3 * i + 2) += -dtmp * rz;
    }
  }
  return EXIT_SUCCESS;
}

}  // namespace dftd
