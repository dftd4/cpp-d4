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
#include <iostream>

#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_multicharge_param.h"
#include "dftd_ncoord.h"

namespace dftd4 {

NCoordBase::NCoordBase(
  double optional_kcn,
  double optional_norm_exp,
  double optional_cutoff,
  double optional_f_directed,
  double optional_cn_max,
  const double *optional_rcov
)
  : kcn(optional_kcn), norm_exp(optional_norm_exp), cutoff(optional_cutoff),
    f_directed(optional_f_directed), cn_max(optional_cn_max),
    rcov(optional_rcov) {}

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
const double NCoordBase::rad[119]{
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

static const double hlfosqrtpi =
  1.0 / 1.7724538509055159; // one over square root of pi

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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

void initializeRealIdx(int nat, TVector<int> &realIdx) {
  realIdx.NewVec(nat);
  for (int i = 0; i < nat; ++i) {
    realIdx(i) = i;
  }
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int NCoordBase::get_ncoord(
  const TMolecule &mol,
  const TMatrix<double> &dist,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  TIVector realIdx;
  initializeRealIdx(mol.NAtoms, realIdx);

  return get_ncoord(mol, realIdx, dist, cn, dcndr, lgrad);
};

int NCoordBase::get_ncoord(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  int info;

  const int nat = realIdx.Max() + 1;
  cn.NewVector(nat);
  if (lgrad) dcndr.NewMatrix(nat, 3 * nat);

  if (lgrad) {
    info = dr_ncoord_base(mol, realIdx, dist, cn, dcndr);
  } else {
    info = ncoord_base(mol, realIdx, dist, cn);
  }
  if (info != EXIT_SUCCESS) return info;

  if (cn_max > 0.0) { // cn_max = -1.0 for EEQ-BC for cn and qloc; using
                      // NCoordErf and NCoordErfEN
    info = cut_coordination_number(cn_max, cn, dcndr, lgrad);
    if (info != EXIT_SUCCESS) return info;
  }

  return EXIT_SUCCESS;
};

int NCoordBase::ncoord_base(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TVector<double> &cn
) {
  double r = 0.0, rcovij = 0.0;
  double countf = 0.0;
  double f_en = 0.0;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      r = dist(ii, jj);
      if (r > cutoff) continue;

      rcovij = get_rcov(mol.ATNO(i), mol.ATNO(j));
      f_en = get_en_factor(mol.ATNO(i), mol.ATNO(j));
      countf = f_en * count_fct(r, rcovij);
      cn(ii) += countf;
      cn(jj) += countf * f_directed;
    }
  }
  return EXIT_SUCCESS;
}

int NCoordBase::dr_ncoord_base(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TVector<double> &cn,
  TMatrix<double> &dcndr
) {
  double r = 0.0, rcovij = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0;
  double countf = 0.0, dcountf = 0.0;
  double f_en = 0;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      r = dist(ii, jj);
      if (r > cutoff) continue;

      rx = (mol.CC(i, 0) - mol.CC(j, 0)) / r;
      ry = (mol.CC(i, 1) - mol.CC(j, 1)) / r;
      rz = (mol.CC(i, 2) - mol.CC(j, 2)) / r;

      rcovij = get_rcov(mol.ATNO(i), mol.ATNO(j));

      f_en = get_en_factor(mol.ATNO(i), mol.ATNO(j));
      countf = f_en * count_fct(r, rcovij);
      cn(ii) += countf;
      cn(jj) += countf * f_directed;

      dcountf = f_en * dr_count_fct(r, rcovij) / pow(rcovij, norm_exp);
      dcndr(jj, 3 * jj) -= dcountf * rx * f_directed;
      dcndr(jj, 3 * jj + 1) -= dcountf * ry * f_directed;
      dcndr(jj, 3 * jj + 2) -= dcountf * rz * f_directed;

      dcndr(jj, 3 * ii) += dcountf * rx * f_directed;
      dcndr(jj, 3 * ii + 1) += dcountf * ry * f_directed;
      dcndr(jj, 3 * ii + 2) += dcountf * rz * f_directed;

      dcndr(ii, 3 * jj) -= dcountf * rx;
      dcndr(ii, 3 * jj + 1) -= dcountf * ry;
      dcndr(ii, 3 * jj + 2) -= dcountf * rz;

      dcndr(ii, 3 * ii) += dcountf * rx;
      dcndr(ii, 3 * ii + 1) += dcountf * ry;
      dcndr(ii, 3 * ii + 2) += dcountf * rz;
    }
  }

  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

double NCoordErf::count_fct(double r, double rc) const {
  return 0.5 * (1.0 + erf(-kcn * (r - rc) / pow(rc, norm_exp)));
}

double NCoordErf::dr_count_fct(double r, double rc) const {
  const double rc_norm_exp = pow(rc, norm_exp);
  const double exponent_term = -pow(kcn * (r - rc) / rc_norm_exp, 2);
  return -kcn * hlfosqrtpi * exp(exponent_term);
}

double NCoordErfEN::count_fct(double r, double rc) const {
  return 0.5 * (1.0 + erf(-kcn * (r - rc) / pow(rc, norm_exp)));
}

double NCoordErfEN::dr_count_fct(double r, double rc) const {
  const double rc_norm_exp = pow(rc, norm_exp);
  const double exponent_term = -pow(kcn * (r - rc) / rc_norm_exp, 2);
  return -kcn * hlfosqrtpi * exp(exponent_term);
}

double NCoordErfD4::count_fct(double r, double rc) const {
  return 0.5 * (1.0 + erf(-kcn * (r - rc) / pow(rc, norm_exp)));
}

double NCoordErfD4::dr_count_fct(double r, double rc) const {
  const double rc_norm_exp = pow(rc, norm_exp);
  const double exponent_term = -pow(kcn * (r - rc) / rc_norm_exp, 2);
  return -kcn * hlfosqrtpi * exp(exponent_term);
}

inline double log_cn_cut(const double cn_max, const double cn) {
  return log(1.0 + exp(cn_max)) - log(1.0 + exp(cn_max - cn));
}

inline double dlog_cn_cut(const double cn_max, const double cn) {
  return exp(cn_max) / (exp(cn_max) + exp(cn));
}

int NCoordErf::cut_coordination_number(
  const double cn_max,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  if (lgrad) {
    // cutting the cn is not (anti)symmetric, so dcndr is not antisymmetric
    // anymore
    double dcnpdcn;
    for (int i = 0; i != cn.N; i++) {
      dcnpdcn = dlog_cn_cut(cn_max, cn(i));
      for (int j = 0; j != cn.N; j++) {
        dcndr(i, 3 * j) *= dcnpdcn;
        dcndr(i, 3 * j + 1) *= dcnpdcn;
        dcndr(i, 3 * j + 2) *= dcnpdcn;
      }
    }
  }

  for (int i = 0; i != cn.N; i++) {
    cn(i) = log_cn_cut(cn_max, cn(i));
  }

  return EXIT_SUCCESS;
}

int NCoordErfEN::cut_coordination_number(
  const double cn_max,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  return EXIT_SUCCESS;
}

int NCoordErfD4::cut_coordination_number(
  const double cn_max,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
) {
  return EXIT_SUCCESS;
}

} // namespace dftd4
