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
 * Definition of the coordination numbers used in DFT-D4 and the
 * electronegativity equilibration (EEQ) model
 *
 * This module works on a distance matrix to avoid recalculating
 * the distances every time.
 */
#include <cmath>
#include <iostream>

#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_ncoord.h"

namespace dftd4 {

/**
 * Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
 * 188-197), values for metals decreased by 10 %.
 *
 * These values are actually never used in the code.
 * Only the scaled values below are used (`rad`).
 */
static const double covalent_rad_d3[119]{
  0.0,  0.32, 0.46,                                           //  H,He
  1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,             //  Li-Ne
  1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,             //  Na-Ar
  1.76, 1.54,                                                 //  K,Ca
  1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09, //  Sc-Zn
  1.12, 1.09, 1.15, 1.10, 1.14, 1.17,                         //  Ga-Kr
  1.89, 1.67,                                                 //  Rb,Sr
  1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, //  Y-Cd
  1.28, 1.26, 1.26, 1.23, 1.32, 1.31,                         //  In-Xe
  2.09, 1.76,                                                 //  Cs,Ba
  1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51,                   //  La-Eu
  1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,                   //  Gd-Yb
  1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32, //  Lu-Hg
  1.30, 1.30, 1.36, 1.31, 1.38, 1.42,                         //  Tl-Rn
  2.01, 1.81,                                                 //  Fr,Ra
  1.67, 1.58, 1.52, 1.53, 1.54, 1.55, 1.49,                   //  Ac-Am
  1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58,                   //  Cm-No
  1.45, 1.41, 1.34, 1.29, 1.27, 1.21, 1.16, 1.15, 1.09, 1.22, //  Lr-Cn
  1.36, 1.43, 1.46, 1.58, 1.48, 1.57                          //  Nh-Og
};

/**
 * D3 covalent radii used to construct the coordination number.
 * rad = covalent_rad_d3 * 4.0/3.0 * aatoau
 *
 * convert bohr (a.u.) to Ångström and back via:
 * static const double autoaa = 0.52917721090449243;
 * static const double aatoau = 1.0 / autoaa;
 */
static double rad[119]{
  0.00000000000000, 0.80628314650472, 1.15903202310054, 3.02356179939270,
  2.36845674285762, 1.94011882127699, 1.88972612462044, 1.78894073130735,
  1.58736994468117, 1.61256629300944, 1.68815533799426, 3.52748876595816,
  3.14954354103407, 2.84718736109480, 2.62042022614034, 2.77159831610998,
  2.57002752948380, 2.49443848449898, 2.41884943951416, 4.43455730577597,
  3.88023764255397, 3.35111432766025, 3.07395449604925, 3.04875814772098,
  2.77159831610998, 2.69600927112516, 2.62042022614034, 2.51963483282725,
  2.49443848449898, 2.54483118115553, 2.74640196778171, 2.82199101276652,
  2.74640196778171, 2.89758005775134, 2.77159831610998, 2.87238370942307,
  2.94797275440789, 4.76210983404351, 4.20779017082151, 3.70386320425606,
  3.50229241762988, 3.32591797933198, 3.12434719270579, 2.89758005775134,
  2.84718736109480, 2.84718736109480, 2.72120561945343, 2.89758005775134,
  3.09915084437752, 3.22513258601889, 3.17473988936234, 3.17473988936234,
  3.09915084437752, 3.32591797933198, 3.30072163100370, 5.26603680060896,
  4.43455730577597, 4.08180842918015, 3.70386320425606, 3.98102303586706,
  3.95582668753879, 3.93063033921052, 3.90543399088224, 3.80464859756915,
  3.82984494589743, 3.80464859756915, 3.77945224924088, 3.75425590091261,
  3.75425590091261, 3.72905955258434, 3.85504129422570, 3.67866685592779,
  3.45189972097334, 3.30072163100370, 3.09915084437752, 2.97316910273616,
  2.92277640607961, 2.79679466443825, 2.82199101276652, 2.84718736109480,
  3.32591797933198, 3.27552528267543, 3.27552528267543, 3.42670337264507,
  3.30072163100370, 3.47709606930161, 3.57788146261470, 5.06446601398278,
  4.56053904741733, 4.20779017082151, 3.98102303586706, 3.82984494589743,
  3.85504129422570, 3.88023764255397, 3.90543399088224, 3.75425590091261,
  3.75425590091261, 3.80464859756915, 3.80464859756915, 3.72905955258434,
  3.77945224924088, 3.93063033921052, 3.98102303586706, 3.65347050759952,
  3.55268511428643, 3.37631067598852, 3.25032893434716, 3.19993623769061,
  3.04875814772098, 2.92277640607961, 2.89758005775134, 2.74640196778171,
  3.07395449604925, 3.42670337264507, 3.60307781094297, 3.67866685592779,
  3.98102303586706, 3.72905955258434, 3.95582668753879,
};

// pauling EN's
static const double en[119]{
  0.0,  // dummy
  2.20, // H
  3.00, // He
  0.98, // Li (2nd)
  1.57, // Be
  2.04, // B
  2.55, // C
  3.04, // N
  3.44, // O
  3.98, // F
  4.50, // Ne
  0.93, // Na (3rd)
  1.31, // Mg
  1.61, // Al
  1.90, // Si
  2.19, // P
  2.58, // S
  3.16, // Cl
  3.50, // Ar
  0.82, // K  (4th)
  1.00, // Ca
  1.36, // Sc
  1.54, // Ti
  1.63, // V
  1.66, // Cr
  1.55, // Mn
  1.83, // Fe
  1.88, // Co
  1.91, // Ni
  1.90, // Cu
  1.65, // Zn
  1.81, // Ga
  2.01, // Ge
  2.18, // As
  2.55, // Se
  2.96, // Br
  3.00, // Kr
  0.82, // Rb (5th)
  0.95, // Sr
  1.22, // Y
  1.33, // Zr
  1.60, // Nb
  2.16, // Mo
  1.90, // Tc
  2.20, // Ru
  2.28, // Rh
  2.20, // Pd
  1.93, // Ag
  1.69, // Cd
  1.78, // In
  1.96, // Sn
  2.05, // Sb
  2.10, // Te
  2.66, // I
  2.60, // Xe
  0.79, // Cs (6th)
  0.89, // Ba
  1.10, // La
  1.12, // Ce
  1.13, // Pr
  1.14, // Nd
  1.15, // Pm
  1.17, // Sm
  1.18, // Eu
  1.20, // Gd
  1.21, // Tb
  1.22, // Dy
  1.23, // Ho
  1.24, // Er
  1.25, // Tm
  1.26, // Yb
  1.27, // Lu
  1.30, // Hf
  1.50, // Ta
  2.36, // W
  1.90, // Re
  2.20, // Os
  2.20, // Ir
  2.28, // Pt
  2.54, // Au
  2.00, // Hg
  1.62, // Tl
  2.33, // Pb
  2.02, // Bi
  2.00, // Po
  2.20, // At
  2.20, // Rn
  0.79, // Fr (7th)
  0.90, // Ra
  1.10, // Ac
  1.30, // Th
  1.50, // Pa
  1.38, // U
  1.36, // Np
  1.28, // Pu
  1.30, // Am
  1.30, // Cm
  1.30, // Bk
  1.30, // Cf
  1.30, // Es
  1.30, // Fm
  1.30, // Md
  1.30, // No
  1.30, // Lr
  1.50, // Rf (only dummies from here)
  1.50, // Db
  1.50, // Sg
  1.50, // Bh
  1.50, // Hs
  1.50, // Mt
  1.50, // Ds
  1.50, // Rg
  1.50, // Cn
  1.50, // Nh
  1.50, // Fl
  1.50, // Lv
  1.50, // Mc
  1.50, // Ts
  1.50, // Og
};

static const double kn = 7.5;
static const double k4 = 4.10451;
static const double k5 = 19.08857;
static const double k6 = 2 * pow(11.28174, 2);
static const double hlfosqrtpi = 1.0 / 1.7724538509055159;

// Maximum CN (not strictly obeyed)
static const double cn_max = 8.0;

int calc_distances(
  const TMolecule &mol,
  const TIVector &realIdx,
  TMatrix<double> &dist
) {
  double rx = 0.0, ry = 0.0, rz = 0.0, tmp = 0.0;
  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;
    dist(ii, ii) = 0.0;

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      rx = mol.CC(i, 0) - mol.CC(j, 0);
      ry = mol.CC(i, 1) - mol.CC(j, 1);
      rz = mol.CC(i, 2) - mol.CC(j, 2);
      tmp = sqrt(rx * rx + ry * ry + rz * rz);
      dist(ii, jj) = tmp;
      dist(jj, ii) = tmp;
    }
  }

  return EXIT_SUCCESS;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int get_ncoord_erf(
  const TMolecule &mol,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  TIVector realIdx;
  initializeRealIdx(mol.NAtoms, realIdx);

  return get_ncoord_erf(mol, realIdx, dist, cutoff, cn, dcndr, lgrad);
};

int get_ncoord_erf(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  int info;

  if (lgrad) {
    info = dncoord_erf(mol, realIdx, dist, cutoff, cn, dcndr);
  } else {
    info = ncoord_erf(mol, realIdx, dist, cutoff, cn);
  }
  if (info != EXIT_SUCCESS) return info;

  info = cut_coordination_number(cn_max, cn, dcndr, lgrad);
  if (info != EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
};

int ncoord_erf(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn
) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double countf = 0.0;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      r = dist(ii, jj);
      if (r > cutoff) continue;

      rcovij = rad[mol.ATNO(i)] + rad[mol.ATNO(j)];
      rr = r / rcovij;
      countf = erf_count(kn, rr);
      cn(ii) += countf;
      cn(jj) += countf;
    }
  }

  return EXIT_SUCCESS;
}

int dncoord_erf(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr
) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0;
  double countf = 0.0, dcountf = 0.0;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      r = dist(ii, jj);
      if (r > cutoff) continue;

      rx = (mol.CC(j, 0) - mol.CC(i, 0)) / r;
      ry = (mol.CC(j, 1) - mol.CC(i, 1)) / r;
      rz = (mol.CC(j, 2) - mol.CC(i, 2)) / r;

      rcovij = rad[mol.ATNO(i)] + rad[mol.ATNO(j)];
      rr = r / rcovij;

      countf = erf_count(kn, rr);
      cn(ii) += countf;
      cn(jj) += countf;

      dcountf = derf_count(kn, rr) / rcovij;
      dcndr(jj, 3 * jj) += dcountf * rx;
      dcndr(jj, 3 * jj + 1) += dcountf * ry;
      dcndr(jj, 3 * jj + 2) += dcountf * rz;

      dcndr(jj, 3 * ii) += dcountf * rx;
      dcndr(jj, 3 * ii + 1) += dcountf * ry;
      dcndr(jj, 3 * ii + 2) += dcountf * rz;

      dcndr(ii, 3 * jj) -= dcountf * rx;
      dcndr(ii, 3 * jj + 1) -= dcountf * ry;
      dcndr(ii, 3 * jj + 2) -= dcountf * rz;

      dcndr(ii, 3 * ii) -= dcountf * rx;
      dcndr(ii, 3 * ii + 1) -= dcountf * ry;
      dcndr(ii, 3 * ii + 2) -= dcountf * rz;
    }
  }

  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int get_ncoord_d4(
  const TMolecule &mol,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  TIVector realIdx;
  initializeRealIdx(mol.NAtoms, realIdx);

  return get_ncoord_d4(mol, realIdx, dist, cutoff, cn, dcndr, lgrad);
};

int get_ncoord_d4(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  if (lgrad) return dncoord_d4(mol, realIdx, dist, cutoff, cn, dcndr);
  return ncoord_d4(mol, realIdx, dist, cutoff, cn);
};

int ncoord_d4(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn
) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double den = 0.0;
  double countf = 0.0;
  int izp, jzp;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    izp = mol.ATNO(i);
    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      r = dist(ii, jj);
      if (r > cutoff) continue;

      jzp = mol.ATNO(j);
      rcovij = rad[izp] + rad[jzp];
      rr = r / rcovij;
      den = k4 * exp(-pow((fabs(en[izp] - en[jzp]) + k5), 2) / k6);
      countf = den * erf_count(kn, rr);

      cn(ii) += countf;
      cn(jj) += countf;
    }
  }
  return EXIT_SUCCESS;
}

int dncoord_d4(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr
) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0;
  double countf = 0.0, dcountf = 0.0, den = 0.0;
  int izp, jzp;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    izp = mol.ATNO(i);
    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      r = dist(ii, jj);
      if (r > cutoff) continue;

      jzp = mol.ATNO(j);
      rx = (mol.CC(j, 0) - mol.CC(i, 0)) / r;
      ry = (mol.CC(j, 1) - mol.CC(i, 1)) / r;
      rz = (mol.CC(j, 2) - mol.CC(i, 2)) / r;

      rcovij = rad[izp] + rad[jzp];
      rr = r / rcovij;
      den = k4 * exp(-pow((fabs(en[izp] - en[jzp]) + k5), 2) / k6);
      countf = den * erf_count(kn, rr);
      cn(ii) += countf;
      cn(jj) += countf;

      dcountf = den * derf_count(kn, rr) / rcovij;
      dcndr(3 * jj, jj) += dcountf * rx;
      dcndr(3 * jj + 1, jj) += dcountf * ry;
      dcndr(3 * jj + 2, jj) += dcountf * rz;
      dcndr(3 * jj, ii) = dcountf * rx;
      dcndr(3 * jj + 1, ii) = dcountf * ry;
      dcndr(3 * jj + 2, ii) = dcountf * rz;
      dcndr(3 * ii, jj) = -dcountf * rx;
      dcndr(3 * ii + 1, jj) = -dcountf * ry;
      dcndr(3 * ii + 2, jj) = -dcountf * rz;
      dcndr(3 * ii, ii) += -dcountf * rx;
      dcndr(3 * ii + 1, ii) += -dcountf * ry;
      dcndr(3 * ii + 2, ii) += -dcountf * rz;
    }
  }
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

double erf_count(double k, double rr) {
  return 0.5 * (1.0 + erf(-k * (rr - 1.0)));
}

double derf_count(double k, double rr) {
  return -k * hlfosqrtpi * exp(-pow(k * (rr - 1.0), 2));
}

int cut_coordination_number(
  const double cn_max,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  if (lgrad) {
    double dcnpdcn;
    for (int i = 0; i != cn.N; i++) {
      dcnpdcn = dlog_cn_cut(cn_max, cn(i));
      for (int j = 0; j != cn.N; j++) {
        dcndr(j, 3 * i) *= dcnpdcn;
        dcndr(j, 3 * i + 1) *= dcnpdcn;
        dcndr(j, 3 * i + 2) *= dcnpdcn;
      }
    }
  }

  for (int i = 0; i != cn.N; i++) {
    cn(i) = log_cn_cut(cn_max, cn(i));
  }

  return EXIT_SUCCESS;
}

inline double log_cn_cut(const double cn_max, const double cn) {
  return log(1.0 + exp(cn_max)) - log(1.0 + exp(cn_max - cn));
}

inline double dlog_cn_cut(const double cn_max, const double cn) {
  return exp(cn_max) / (exp(cn_max) + exp(cn));
}

} // namespace dftd4
