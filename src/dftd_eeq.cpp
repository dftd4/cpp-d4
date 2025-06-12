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
 * Electronegativity equilibration (EEQ) model for DFT-D4.
 * This implementation contains only the essential parts for DFT-D4.
 */
#include <cmath>

#include "dftd_cblas.h"
#include "dftd_eeq.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_ncoord.h"
#include "dftd_eeq_param.h"

// wrap all charge models in the multicharge namespace
namespace multicharge {

// constants
static const double small = 1e-14;
static const double pi = 3.1415926535897932384626433832795029;
static const double sqrtpi = std::sqrt(pi);
static const double sqrt2pi = std::sqrt(2.0 / pi);

// Base class for charge models
ChargeModel::ChargeModel(){}

// Get charges after adjusting atom indices in case ghost atoms are present
int ChargeModel::get_charges(
  const TMolecule &mol,
  const TMatrix<double> &dist,
  const int charge,
  const double cutoff,
  TVector<double> &q,
  TMatrix<double> &dqdr,
  bool lgrad
) {
  TIVector realIdx;
  initializeRealIdx(mol.NAtoms, realIdx);

  return get_charges(mol, realIdx, dist, charge, cutoff, q, dqdr, lgrad);
};

// Get charges
int ChargeModel::get_charges(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const int charge,
  const double cutoff,
  TVector<double> &q,
  TMatrix<double> &dqdr,
  bool lgrad
) {
  int info{0};
  bool lverbose{false};
  int nat = realIdx.Max() + 1;

  TVector<double> cn;    // EEQ cordination number
  TMatrix<double> dcndr; // Derivative of EEQ-CN

  cn.NewVec(nat);
  if (lgrad) dcndr.NewMat(nat, 3 * nat);

  // get the EEQ coordination number
  info = get_ncoord_erf(mol, realIdx, dist, cutoff, cn, dcndr, lgrad);
  if (info != EXIT_SUCCESS) return info;

  // corresponds to model%solve in Fortran
  info =
    eeq_chrgeq(mol, realIdx, dist, charge, cn, q, dcndr, dqdr, lgrad, lverbose);
  if (info != EXIT_SUCCESS) return info;

  dcndr.DelMat();
  cn.DelVec();

  return EXIT_SUCCESS;
};

int ChargeModel::eeq_chrgeq(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const int &charge,
  const TVector<double> &cn,
  TVector<double> &qvec,
  TMatrix<double> &dcndr,
  TMatrix<double> &dqdr,
  bool lgrad /*= false*/,
  bool lverbose /*= false*/
) {
  double qtotal = 0.0;
  int info = 0;
  int n = realIdx.Max() + 1;
  int m = n + 1;

  TMatrix<double> Amat; // Coulomb matrix
  TVector<double> xvec; // x (chi) vector
  Amat.NewMat(m, m);
  xvec.NewVec(m);

  TVector<double> dxdcn; // Derivative of chi vector w.r.t. CN
  if (lgrad) dxdcn.NewVec(m);

  info = get_vrhs(mol, realIdx, charge, cn, xvec, dxdcn, lgrad);
  if (info != EXIT_SUCCESS) return info;

  info = get_amat_0d(mol, realIdx, dist, Amat);
  if (info != EXIT_SUCCESS) return info;

  TVector<double> vrhs;
  vrhs.NewVec(m);

  TMatrix<double> Ainv;
  Ainv.NewMat(m, m);
  Ainv.CopyMat(Amat);

  // solve: A Q = X (Eq.4) -> Q = Ainv X
  info = BLAS_InvertMatrix(Ainv);
  if (info != EXIT_SUCCESS) return info;

  // no return in ORCA
  BLAS_Add_Mat_x_Vec(vrhs, Ainv, xvec, false, 1.0);
  // if (info != EXIT_SUCCESS) return info;

  // remove charge constraint (make vector smaller by one) and total charge
  qtotal = 0.0;
  for (int iat = 0, ii = 0; iat != mol.NAtoms; iat++) {
    ii = realIdx(iat);
    if (ii < 0) continue;

    qvec(ii) = vrhs(ii);
    qtotal += qvec(ii);
  }

  // check total charge and additional printout
  if (fabs(qtotal - charge) > 1.0e-8) {
    printf(
      "DFT-D4: EEQ charge constraint error: %14.8f vs. %14d\n", qtotal, charge
    );
  }

  if (lverbose) {
    printf("    #   sym             EN              q            Aii\n");
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;
      printf(
        "%5d %5d %14.8f %14.8f %14.8f\n",
        i,
        mol.ATNO(i),
        -xvec(ii),
        qvec(ii),
        Amat(ii, ii)
      );
    }
  }

  // Gradient (note that the corresponding gradient flag in Fortran is `cpq`)
  if (lgrad) {
    TMatrix<double> dAmat;
    dAmat.NewMat(3 * n, m);
    TMatrix<double> atrace;
    atrace.NewMat(m, 3);

    info = get_damat_0d(mol, realIdx, dist, vrhs, Amat, dAmat, atrace);
    if (info != EXIT_SUCCESS) return info;

    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      dAmat(3 * ii, ii) += atrace(ii, 0);
      dAmat(3 * ii + 1, ii) += atrace(ii, 1);
      dAmat(3 * ii + 2, ii) += atrace(ii, 2);

      for (int j = 0, jj = 0; j != mol.NAtoms; j++) {
        jj = realIdx(j);
        if (jj < 0) continue;

        dAmat(3 * jj, ii) -= dcndr(jj, 3 * ii) * dxdcn(ii);
        dAmat(3 * jj + 1, ii) -= dcndr(jj, 3 * ii + 1) * dxdcn(ii);
        dAmat(3 * jj + 2, ii) -= dcndr(jj, 3 * ii + 2) * dxdcn(ii);
      }
    }

    // we do not need these gradient-related matrices anymore
    atrace.DelMat();
    dxdcn.DelVec();

    // Ainv with last column removed
    TMatrix<double> A;
    A.NewMat(Ainv.rows, Ainv.cols - 1);
    for (int i = 0; i < Ainv.rows; i++) {
      for (int j = 0; j < Ainv.cols - 1; j++) {
        A(i, j) = Ainv(i, j);
      }
    }

    // no return in ORCA
    BLAS_Add_Mat_x_Mat(dqdr, dAmat, A, false, false, -1.0);
    // if (info != EXIT_SUCCESS) return info;

    dAmat.DelMat();
  }

  // free all memory
  Ainv.DelMat();
  Amat.DelMat();
  xvec.DelVec();

  return EXIT_SUCCESS;
}

// EEQ Model class derived from ChargeModel
EEQModel::EEQModel()
  // Constructor initializer list
  : xi(EEQ_param::xi),
    gam(EEQ_param::gam),
    kappa(EEQ_param::kappa),
    alp(EEQ_param::alp)
{}

int EEQModel::get_vrhs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int &charge,
  const TVector<double> &cn,
  TVector<double> &Xvec,
  TVector<double> &dXvec,
  bool lgrad
) const {
  double tmp{0.0};
  int izp;
  int nat = realIdx.Max() + 1;

  if (lgrad) {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(cn(ii) + small);
      Xvec(ii) = -xi[izp] + tmp * cn(ii);
      dXvec(ii) = 0.5 * tmp;
    }
    dXvec(nat) = 0.0;
  } else {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(cn(ii) + small);
      Xvec(ii) = -xi[izp] + tmp * cn(ii);
    }
  }

  // place charge at last index of xvec
  Xvec(nat) = charge;

  return EXIT_SUCCESS;
};

int EEQModel::get_amat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TMatrix<double> &Amat
) const {
  double gamij = 0.0;
  int mm = realIdx.Max() + 1;
  int izp, jzp;
  double alphai, alphaj;
  double tmp, r;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    izp = mol.ATNO(i);
    alphai = pow(alp[izp], 2);
    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      jzp = mol.ATNO(j);
      alphaj = pow(alp[jzp], 2);

      r = dist(ii, jj);
      gamij = 1.0 / std::sqrt(alphai + alphaj);
      tmp = std::erf(gamij * r) / r;
      Amat(ii, jj) = tmp;
      Amat(jj, ii) = tmp;
    }
    gamij = gam[izp];
    Amat(ii, ii) = gamij + sqrt2pi / alp[izp];
    Amat(ii, mm) = 1.0;
    Amat(mm, ii) = 1.0;
  }
  Amat(mm, mm) = 0.0;

  return EXIT_SUCCESS;
};

int EEQModel::get_damat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const TVector<double> &q,
  const TMatrix<double> &Amat,
  TMatrix<double> &dAmat,
  TMatrix<double> &atrace
) const {
  double alphai, alphaj;
  double rx, ry, rz, r2;
  double arg, gam, dtmp;
  double dgx, dgy, dgz;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    alphai = pow(alp[mol.ATNO(i)], 2);

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      alphaj = pow(alp[mol.ATNO(j)], 2);

      rx = mol.CC(i, 0) - mol.CC(j, 0);
      ry = mol.CC(i, 1) - mol.CC(j, 1);
      rz = mol.CC(i, 2) - mol.CC(j, 2);
      r2 = pow(dist(ii, jj), 2);

      gam = 1.0 / std::sqrt((alphai + alphaj));
      arg = gam * gam * r2;
      dtmp = 2.0 * gam * std::exp(-arg) / (sqrtpi * r2) - Amat(jj, ii) / r2;
      dgx = dtmp * rx;
      dgy = dtmp * ry;
      dgz = dtmp * rz;

      atrace(ii, 0) += dgx * q(jj);
      atrace(ii, 1) += dgy * q(jj);
      atrace(ii, 2) += dgz * q(jj);
      atrace(jj, 0) -= dgx * q(ii);
      atrace(jj, 1) -= dgy * q(ii);
      atrace(jj, 2) -= dgz * q(ii);

      dAmat(3 * ii, jj) = dgx * q(ii);
      dAmat(3 * ii + 1, jj) = dgy * q(ii);
      dAmat(3 * ii + 2, jj) = dgz * q(ii);
      dAmat(3 * jj, ii) = -dgx * q(jj);
      dAmat(3 * jj + 1, ii) = -dgy * q(jj);
      dAmat(3 * jj + 2, ii) = -dgz * q(jj);
    }
  }

  return EXIT_SUCCESS;
};


} // namespace multicharge
