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
#include "dftd_multicharge_param.h"

// wrap all charge models in the multicharge namespace
namespace multicharge {
  using dftd4::NCoordErf;
  using dftd4::NCoordErfEN;

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
  dftd4::initializeRealIdx(mol.NAtoms, realIdx);

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
  const int nat = realIdx.Max() + 1;
  TVector<double> cn;  // coordination number
  TMatrix<double> dcndr;  // derivative of the coordination number

  // get correct cn for current charge model
  get_cn(mol, realIdx, dist, cn, dcndr, lgrad);

  // corresponds to model%solve in Fortran
  info =
    eeq_chrgeq(mol, realIdx, dist, cn, dcndr, charge, q, dqdr, lgrad, lverbose);
  if (info != EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
};

int ChargeModel::eeq_chrgeq(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const TVector<double> &cn,
  const TMatrix<double> &dcndr,
  const int &charge,
  TVector<double> &qvec,
  TMatrix<double> &dqdr,
  bool lgrad /*= false*/,
  bool lverbose /*= false*/
) {
  int info{0};
  const int n = realIdx.Max() + 1;  // Number of atoms
  int m = n + 1;  // Number of atoms plus one for constraint

  TMatrix<double> Amat;  // Coulomb matrix
  TVector<double> xvec;  // x (chi) vector, practically the right-hand side of the set of linear equations
  TMatrix<double> dxvecdr;  // Derivative of the x vector w.r.t. atom positions
  Amat.NewMat(m, m);
  xvec.NewVec(m);

  // EEQ-BC specific variables:
  TVector<double> qloc;     // Local charge (allocated in get_qloc)
  TMatrix<double> dqlocdr;  // Derivative of local charge w.r.t. atom positions (allocated in get_ncoord)
  TMatrix<double> cmat;     // Capacitance matrix (allocated in get_cmat)
  TMatrix<double> dcmatdr;  // Derivative of the capacitance matrix w.r.t. atom positions (allocated in get_dcmatdr)

  info = get_vrhs(mol, realIdx, charge, dist, cn, dcndr, xvec, dxvecdr, qloc, dqlocdr, cmat, dcmatdr, lgrad);
  if (info != EXIT_SUCCESS) return info;

  info = get_amat_0d(mol, realIdx, dist, cn, qloc, cmat, Amat);
  if (info != EXIT_SUCCESS) return info;

  TVector<double> vrhs;  // right-hand side (RHS) of linear equations; finally holds the solution of the linear equations
  vrhs.NewVec(m);

  TMatrix<double> Ainv;  // Inverse of the Coulomb matrix for dqdr
  // Solve the set of linear equations with the Inverse Ainv
  //   which we need later for dqdr as requested by lgrad=true
  if (lgrad) {
    Ainv.NewMat(m, m);
    Ainv.CopyMat(Amat);

    // solve: A Q = X
    info = BLAS_InvertMatrix(Ainv);
    if (info != EXIT_SUCCESS) return info;

    // no return in ORCA
    BLAS_Add_Mat_x_Vec(vrhs, Ainv, xvec, false, 1.0);
    // if (info != EXIT_SUCCESS) return info;
  }
  // Solve the set of linear equations directly (faster since we do not need Ainv)
  else {
    TVector<double> rhs;
    rhs.NewVec(m);
    rhs.CopyVec(xvec);

    // Now solve (A * q = x) for symmetric Coulomb matrix Amat
    // Note that only the lower triangle of Amat is updated
    info = BLAS_SolveSymmetric(Amat, rhs);
    if (info != EXIT_SUCCESS) return info;

    // Extract the solution
    for (int i = 0, ii = 0; i < mol.NAtoms; i++) {
        ii = realIdx(i);
        if (ii < 0) continue;
        vrhs(ii) = rhs(ii);
    }
  }

  // Remove charge constraint (make vector smaller by one) and total charge
  double qtotal = 0.0;  // total charge for checking charge constraint
  for (int iat = 0, ii = 0; iat != mol.NAtoms; iat++) {
    ii = realIdx(iat);
    if (ii < 0) continue;

    qvec(ii) = vrhs(ii);  // Assign partial charges to the correct variable
    qtotal += qvec(ii);  // Total charge of the system
  }

  // Check total charge and additional printout
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
    TMatrix<double> dAmatdr;
    dAmatdr.NewMat(3 * n, m);
    TMatrix<double> atrace;
    atrace.NewMat(m, 3);

    // Calculate the derivative of the Coulomb matrix w.r.t. atom positions (dAmatdr)
    info = get_damat_0d(mol, realIdx, dist, vrhs, Amat, dAmatdr, atrace,
       cn, dcndr, qloc, dqlocdr, cmat, dcmatdr);
    if (info != EXIT_SUCCESS) return info;

    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      dAmatdr(3 * ii, ii) += atrace(ii, 0);
      dAmatdr(3 * ii + 1, ii) += atrace(ii, 1);
      dAmatdr(3 * ii + 2, ii) += atrace(ii, 2);

      for (int j = 0, jj = 0; j != mol.NAtoms; j++) {
        jj = realIdx(j);
        if (jj < 0) continue;

        // dxvecdr = del C/del R * X + del X / del R * C 
        // CN contribution
        dAmatdr(3 * jj,     ii) -= dxvecdr(ii, 3 * jj    );
        dAmatdr(3 * jj + 1, ii) -= dxvecdr(ii, 3 * jj + 1);
        dAmatdr(3 * jj + 2, ii) -= dxvecdr(ii, 3 * jj + 2);
      }
    }

    // we do not need these gradient-related matrices anymore
    atrace.DelMat();

    // Ainv with last column removed
    TMatrix<double> A;
    A.NewMat(Ainv.rows, Ainv.cols - 1);
    for (int i = 0; i < Ainv.rows; i++) {
      for (int j = 0; j < Ainv.cols - 1; j++) {
        A(i, j) = Ainv(i, j);
      }
    }

    // Calculate the charge gradient w.r.t. atom positions (dqdr)
    // no return in ORCA
    BLAS_Add_Mat_x_Mat(dqdr, dAmatdr, A, false, false, -1.0);
    // if (info != EXIT_SUCCESS) return info;

    dAmatdr.DelMat();
    Ainv.DelMat();
    A.DelMat();
    dqlocdr.DelMat();
    dxvecdr.DelMat();
    dcmatdr.DelMat();
  }

  // free all memory
  Amat.DelMat();
  xvec.DelVec();
  qloc.DelVec();
  cmat.DelMat();
  vrhs.DelVec();

  return EXIT_SUCCESS;
}

// EEQ Model class derived from ChargeModel
EEQModel::EEQModel()
  // Constructor initializer list
  : xi(multicharge_param::eeq::xi),
    gam(multicharge_param::eeq::gam),
    kappa(multicharge_param::eeq::kappa),
    alp(multicharge_param::eeq::alp)
{
  ncoord_erf = NCoordErf();
}

int EEQModel::get_vrhs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int &charge,
  const TMatrix<double> &dist,
  const TVector<double> &cn,
  const TMatrix<double> &dcndr,
  TVector<double> &xvec,
  TMatrix<double> &dxvecdr,
  TVector<double> &qloc,
  TMatrix<double> &dqlocdr,
  TMatrix<double> &cmat,
  TMatrix<double> &dcmatdr,
  bool lgrad
) {
  TVector<double> dxdcn;
  double tmp{0.0};
  int izp;
  const int nat = realIdx.Max() + 1;

  if (lgrad) {
    dxvecdr.NewMat(nat, 3*nat);
    dxdcn.NewVec(nat);
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(cn(ii) + small);
      xvec(ii) = -xi[izp] + tmp * cn(ii);
      dxdcn(ii) = 0.5 * tmp;

      for (int j = 0, jj = 0; j != mol.NAtoms; j++) {
        jj = realIdx(j);
        if (jj < 0) continue;
        dxvecdr(ii, 3*jj  ) = dxdcn(ii) * dcndr(ii, 3 * jj    ); 
        dxvecdr(ii, 3*jj+1) = dxdcn(ii) * dcndr(ii, 3 * jj + 1); 
        dxvecdr(ii, 3*jj+2) = dxdcn(ii) * dcndr(ii, 3 * jj + 2); 
      }
    }
  } else {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(cn(ii) + small);
      xvec(ii) = -xi[izp] + tmp * cn(ii);
    }
  }

  // place charge at last index of xvec
  xvec(nat) = charge;

  return EXIT_SUCCESS;
};

int EEQModel::get_amat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const TVector<double> &cn,
  const TVector<double> &qloc,
  const TMatrix<double> &cmat,
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
    alphai = alp[izp] * alp[izp];
    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      jzp = mol.ATNO(j);
      alphaj = alp[jzp] * alp[jzp];

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
  TMatrix<double> &dAmatdr,
  TMatrix<double> &atrace,
  const TVector<double> &cn,
  const TMatrix<double> &dcndr,
  const TVector<double> &qloc,
  const TMatrix<double> &dqlocdr,
  const TMatrix<double> &cmat,
  const TMatrix<double> &dcmatdr
) const {
  double alphai, alphaj;
  double rx, ry, rz, r2;
  double arg, gam, dtmp;
  double dgx, dgy, dgz;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    alphai = alp[mol.ATNO(i)] * alp[mol.ATNO(i)];

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      alphaj = alp[mol.ATNO(j)] * alp[mol.ATNO(j)];

      rx = mol.CC(i, 0) - mol.CC(j, 0);
      ry = mol.CC(i, 1) - mol.CC(j, 1);
      rz = mol.CC(i, 2) - mol.CC(j, 2);
      r2 = dist(ii, jj) * dist(ii, jj);

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

      dAmatdr(3 * ii,     jj) =  dgx * q(ii);
      dAmatdr(3 * ii + 1, jj) =  dgy * q(ii);
      dAmatdr(3 * ii + 2, jj) =  dgz * q(ii);
      dAmatdr(3 * jj,     ii) = -dgx * q(jj);
      dAmatdr(3 * jj + 1, ii) = -dgy * q(jj);
      dAmatdr(3 * jj + 2, ii) = -dgz * q(jj);
    }
  }

  return EXIT_SUCCESS;
};

// Calculate the coordination number, forwarding to get_ncoord
int EEQModel::get_cn(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
)
{
  int info{0};

  // get the EEQ coordination number
  info = ncoord_erf.get_ncoord(mol, realIdx, dist, cn, dcndr, lgrad);
  if (info != EXIT_SUCCESS) {
    throw std::runtime_error("EEQModel::get_cn: Failed to compute coordination numbers.");
  }

  return info;
};

// EEQ-BC model class derived from ChargeModel
EEQBCModel::EEQBCModel()
  // Constructor initializer list
  : kcnrad(multicharge_param::eeqbc::kcnrad),
    kbc(multicharge_param::eeqbc::kbc),
    cutoff(multicharge_param::eeqbc::cutoff),
    cn_exp(multicharge_param::eeqbc::cn_exp),
    norm_exp(multicharge_param::eeqbc::norm_exp),
    chi(multicharge_param::eeqbc::eeqbc_chi),
    eta(multicharge_param::eeqbc::eeqbc_eta),
    rad(multicharge_param::eeqbc::eeqbc_rad),
    kcnchi(multicharge_param::eeqbc::eeqbc_kcnchi),
    kqchi(multicharge_param::eeqbc::eeqbc_kqchi),
    kqeta(multicharge_param::eeqbc::eeqbc_kqeta),
    cap(multicharge_param::eeqbc::eeqbc_cap),
    cov_radii(multicharge_param::eeqbc::eeqbc_cov_radii),
    avg_cn(multicharge_param::eeqbc::eeqbc_avg_cn),
    rvdw(multicharge_param::eeqbc::eeqbc_rvdw)
{
  ncoord_erf = NCoordErf(
    ncoorderf_kcn,
    ncoorderf_norm_exp,
    ncoorderf_cutoff,
    ncoorderf_f_directed,
    ncoorderf_cn_max,
    multicharge_param::eeqbc::eeqbc_cov_radii);
}

// Calculate the right-hand side (RHS) vector of the system of linear equations for EEQ-BC
int EEQBCModel::get_vrhs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int &charge,
  const TMatrix<double> &dist,
  const TVector<double> &cn,
  const TMatrix<double> &dcndr,
  TVector<double> &xvec,
  TMatrix<double> &dxvecdr,
  TVector<double> &qloc,
  TMatrix<double> &dqlocdr,
  TMatrix<double> &cmat,
  TMatrix<double> &dcmatdr,
  bool lgrad
) {
  int info{0};
  const int n_atoms = realIdx.Max() + 1;
  // calculate the capacitance matrix
  info = get_cmat(mol, realIdx, dist, cmat);
  if (info != EXIT_SUCCESS) {
     printf("EEQBCModel::get_vrhs: Failed to calculate the capacitance matrix.");
     return info;
  }
  
  if ( ! lgrad) {
    // calculate the right-hand side (RHS)
    info = get_xvec(mol, realIdx, dist, cn, cmat, charge, qloc, xvec);
  } else {
    // calculate RHS and derivative
    info = get_xvec_derivs(mol, realIdx, dist, cn, dcndr, cmat, charge, xvec, dxvecdr, qloc, dqlocdr, dcmatdr);
  }
  if (info != EXIT_SUCCESS) {
     printf("EEQBCModel::get_vrhs: Failed to calculate the right hand side.");
     return info;
  }

  return EXIT_SUCCESS;
};

// Calculate the Coulomb matrix for EEQ-BC
int EEQBCModel::get_amat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const TVector<double> &cn,
  const TVector<double> &qloc,
  const TMatrix<double> &cmat,
  TMatrix<double> &Amat
) const {
  const int nat = realIdx.Max() + 1;
  int iat, jat; // atomic numbers
  double norm_cn; // coordination number normalization factor
  double r, radi, radj, gamij2, tmp;


  for (int i = 0, ii = 0; i < mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;
    iat = mol.ATNO(i);
    norm_cn = 1.0 / pow(avg_cn[iat], norm_exp);
    radi = rad[iat] * (1.0 - kcnrad*cn(ii)*norm_cn);

    for (int j = 0, jj = 0; j < i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;
      jat = mol.ATNO(j);
      r = dist(ii, jj);
      norm_cn =  1.0 / pow(avg_cn[jat], norm_exp);
      radj = rad[jat] * (1.0 - kcnrad*cn(jj)*norm_cn);
      gamij2 = 1.0 / (radi * radi + radj * radj);
      tmp = cmat(jj, ii) * erf(r * sqrt(gamij2)) / r;
      Amat(jj, ii) = Amat(jj, ii) + tmp;
      Amat(ii, jj) = Amat(ii, jj) + tmp;
    }
    tmp = eta[iat] + kqeta[iat] * qloc(ii) + sqrt(2.0/pi) / radi;
    Amat(ii, ii) = Amat(ii, ii) + tmp * cmat(ii, ii) + 1.0;
    // Set entries for charge constraint
    Amat(nat, ii) = 1.0;
    Amat(ii, nat) = 1.0;
  }
  Amat(nat, nat) = 0.0; 

  return EXIT_SUCCESS;
};

// Calculate the derivative of the Coulomb matrix w.r.t. the coordinates for EEQ-BC
int EEQBCModel::get_damat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const TVector<double> &q,
  const TMatrix<double> &Amat,
  TMatrix<double> &dAmatdr,
  TMatrix<double> &atrace,
  const TVector<double> &cn,
  const TMatrix<double> &dcndr,
  const TVector<double> &qloc,
  const TMatrix<double> &dqlocdr,
  const TMatrix<double> &cmat,
  const TMatrix<double> &dcmatdr
) const {

  int izp, jzp;
  double norm_cn;
  double gam;
  double arg, dtmp1, dtmp2, dtmp3, dtmp4, dtmp5;
  double radi, dradcn_i, radj, dradcn_j, r2;
  const int n_atoms = realIdx.Max() + 1;
  TVector<double> dgamdr;
  TVector<double> vec, dG;
  vec.NewVec(3);
  dG.NewVec(3);
  dgamdr.NewVec(3*n_atoms);


  for (int i = 0, ii = 0; i < mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;
    izp = mol.ATNO(i);
    // Effective charge width of i
    norm_cn = 1.0 / pow(avg_cn[izp], norm_exp);
    radi = rad[izp] * (1.0 - kcnrad * cn(ii) * norm_cn);
    dradcn_i = - rad[izp] * kcnrad * norm_cn;
    for (int j = 0, jj = 0; j < i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;
      // Effective charge width of j
      jzp = mol.ATNO(j);
      norm_cn = 1.0 / pow(avg_cn[jzp], norm_exp);
      radj = rad[jzp] * (1.0 - kcnrad * cn(jj) * norm_cn);
      dradcn_j = - rad[jzp] * kcnrad * norm_cn;
      gam = 1.0 / (pow(radi,2) + pow(radj,2));
      for (int k = 0, kk = 0; k < mol.NAtoms; k++) {
        kk = realIdx(k);
        if (kk < 0) continue;
        // Coulomb interaction of Gaussian charges
        for (int c = 0; c < 3; c++) {
          vec(c) = mol.CC(jj, c) - mol.CC(ii, c);
          dgamdr(3*kk+c) =  -(radi * dradcn_i * dcndr(ii, 3*kk + c) + radj * dradcn_j * dcndr(jj, 3*kk + c)) * pow(gam, 3);
        }
      }
      // Explicit derivative
      r2 = pow(vec(0), 2) + pow(vec(1), 2) + pow(vec(2), 2);
      arg = gam * gam * r2;
      dtmp1 = 2.0 * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 *sqrt(r2)); 
      // Effective charge width derivative
      dtmp2 = 2.0 * exp(-arg) / sqrtpi;
      // Capacitance derivative off-diagonal
      dtmp3 = erf(sqrt(r2) * gam) / sqrt(r2);
      // Capacitance derivative diagonal
      dtmp4 = (eta[izp] + kqeta[izp] * qloc(ii) + sqrt2pi / radi) * q(ii);
      dtmp5 = (eta[jzp] + kqeta[jzp] * qloc(jj) + sqrt2pi / radj) * q(jj);
      for (int c = 0; c < 3; c++) {
        // Explicit derivative
        dG(c) = dtmp1 * vec(c);
        atrace(ii, c) += -dG(c) * q(jj) * cmat(ii, jj);
        atrace(jj, c) += +dG(c) * q(ii) * cmat(jj, ii);
        dAmatdr(3*ii+c, jj) += -dG(c) * q(ii) * cmat(jj, ii);
        dAmatdr(3*jj+c, ii) += +dG(c) * q(jj) * cmat(ii, jj);
        // Effective charge width derivative
        atrace(ii, c) += -dtmp2 * q(jj) * dgamdr(3*jj+c) * cmat(ii, jj);
        atrace(jj, c) += -dtmp2 * q(ii) * dgamdr(3*ii+c) * cmat(jj, ii);
        dAmatdr(3*ii+c, jj) += +dtmp2 * q(ii) * dgamdr(3*ii+c) * cmat(jj, ii);
        dAmatdr(3*jj+c, ii) += +dtmp2 * q(jj) * dgamdr(3*jj+c) * cmat(ii, jj);
        // Capacitance derivative off-diagonal
        atrace(ii, c) += -dtmp3 * q(jj) * dcmatdr(ii, 3*jj+c);
        atrace(jj, c) += -dtmp3 * q(ii) * dcmatdr(jj, 3*ii+c);
        dAmatdr(3*ii+c, jj) += +dtmp3 * q(ii) * dcmatdr(jj, 3*ii+c);
        dAmatdr(3*jj+c, ii) += +dtmp3 * q(jj) * dcmatdr(ii, 3*jj+c);
        // Capacitance derivative diagonal
        dAmatdr(3*jj+c, ii) += -dtmp4 * dcmatdr(ii, 3*jj+c);
        dAmatdr(3*ii+c, jj) += -dtmp5 * dcmatdr(jj, 3*ii+c);
      }
    }  // jj
      dtmp1 = kqeta[izp] * q(ii) * cmat(ii, ii);  // Hardness
      dtmp2 = -sqrt2pi * dradcn_i / pow(radi, 2) * q(ii) * cmat(ii, ii);  // Effective charge width
      for (int k = 0, kk = 0; k < mol.NAtoms; k++) {
        kk = realIdx(k);
        if (kk < 0) continue;
        for (int c = 0; c < 3; c++) {
          // Hardness derivative
          dAmatdr(3*kk+c, ii) += +dtmp1 * dqlocdr(ii, 3*kk+c);
          // Effective charge width derivative
          dAmatdr(3*kk+c, ii) += +dtmp2 * dcndr(ii, 3*kk+c);
        }
      }
      // Capacitance derivative
      dtmp3 = (eta[izp] + kqeta[izp] * qloc(ii) + sqrt2pi / radi) * q(ii);
      for (int c = 0; c < 3; c++) {
        dAmatdr(3*ii+c, ii) += +dtmp3 * dcmatdr(ii, 3*ii+c);
      }
}  // ii

  
  return EXIT_SUCCESS;
};

// Get purely geometry-dependent local charges
int EEQBCModel::get_qloc(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const double q_tot,  // total system charge
  TVector<double> &qloc,  // Local charge
  TMatrix<double> &dqlocdr,  // Derivative of local charge w.r.t. atom positions
  bool lgrad
) 
{
  const double cutoff = 25.0;
  TVector<double> cn;
  // Electronegativity scaled coordination number with EEQ-BC parameters
  NCoordErfEN ncoord_erf_en(
    ncoorderf_kcn,
    ncoorderf_norm_exp,
    ncoorderf_cutoff,
    ncoorderfen_f_directed,
    ncoorderf_cn_max);
  const int n_atoms = realIdx.Max() + 1;  // Number of atoms
  qloc.NewVector(n_atoms);
  const double q_tot_norm = q_tot/n_atoms;

  ncoord_erf_en.get_ncoord(mol, dist, cn, dqlocdr, lgrad);

  for (int i = 0, ii = 0; i < mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;
    qloc(ii) = cn(ii) + q_tot_norm;
  }

  return EXIT_SUCCESS;
};


// Calculate the coordination number, forwarding to get_ncoord
int EEQBCModel::get_cn(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
)
{
  int info{0};

  info = ncoord_erf.get_ncoord(mol, realIdx, dist, cn, dcndr, lgrad);
  if (info != EXIT_SUCCESS) {
    throw std::runtime_error("EEQBCModel::get_cn: Failed to compute coordination numbers.");
  }

  return info;
};

// Get capacitance for bond between atoms i and j for EEQ-BC 
int EEQBCModel::get_cpair(
  int iat,  // atom type of i
  int jat,  // atom type of j
  double &dist_ij,  // distance between atom j to atom i
  double &c_ij  // Out: Capacitance for bond ij
) const {
  int iat_zero = iat;  // convert to index counting from zero (e.g. H=0)
  int jat_zero = jat;  // convert to index counting from zero (e.g. H=0)
  int ij_at;
  if (iat_zero > jat_zero) {
    ij_at = (jat_zero + iat_zero*(iat_zero - 1)/2) - 1;
  } else {
    ij_at = (iat_zero + jat_zero*(jat_zero - 1)/2) - 1;
  }
  // Calculate the argument of the error function
  double arg = - kbc*(dist_ij - rvdw[ij_at])/(rvdw[ij_at]);
  c_ij = sqrt(cap[iat]*cap[jat])*0.5*(1.0 + erf(arg));

  return EXIT_SUCCESS;
}

// Get derivative of the capacitance for bond between atoms i and j for EEQ-BC 
int EEQBCModel::get_dcpair(
  double dist_ij,  // distance between atom j to atom i
  double rvdw_ijat,  // pairwise van der Waals radii for ij atom types
  double cap_ij,  // product of bond capacitances for atom types of i and j
  TVector<double> &vec, // Vector from j to i
  TVector<double> &dcdr_ij  // Out: Capacitance for bond ij
) const {
  // Calculate the argument of the error function
  //     arg = -(kbc * (r1      - rvdw       ) / rvdw)**2
  double arg =   kbc * (dist_ij - rvdw_ijat) / rvdw_ijat;
  //     dtmp = sqrt(capi * capj) * kbc * exp(arg) / (sqrtpi * rvdw)
  double dtmp = sqrt(cap_ij) * kbc * exp( - pow(arg, 2)) / (sqrt(pi) * rvdw_ijat);
  for (int c = 0; c < 3; c++) {
    dcdr_ij(c) = dtmp * vec(c) / dist_ij;
  }

  return EXIT_SUCCESS;
}

// Get the capacitance matrix
int EEQBCModel::get_cmat(
  const TMolecule &mol,  // molecular geometry
  const TIVector &realIdx,  // The real atom indices (for excluding dummy atoms)
  const TMatrix<double> &dist, // atom distances
  TMatrix<double> &cmat  // Out: Capacitance matrix
) {
  int iat;  // atom type of i
  int jat;  // atom type of j
  double c_ij;  // Capacitance for bond between atoms i and j
  double dist_ij; // distance between atoms i and j
  const int n_atoms = realIdx.Max() + 1;
  cmat.NewMatrix(n_atoms + 1, n_atoms + 1);
  for (int i = 0, ii = 0; i < mol.NAtoms; i++)
  {
    ii = realIdx(i);
    if (ii < 0) continue;
    iat = mol.ATNO(i);
    for (int j = 0, jj = 0; j < i; j++)
    {
      jj = realIdx(j);
      if (jj < 0) continue;
      jat = mol.ATNO(j);
      dist_ij = dist(ii,jj);
      get_cpair(iat, jat, dist_ij, c_ij);
      // Calulate Off-diagonal elements; bond capacitances
      cmat(ii, jj) = - c_ij;
      cmat(jj, ii) = - c_ij;
      // Calculate diagonal elements; self-capacitance as the negative sum of bond capacitances
      cmat(ii, ii) = cmat(ii, ii) + c_ij;
      cmat(jj, jj) = cmat(jj, jj) + c_ij;
    }
  }
  cmat(n_atoms, n_atoms) = 1.0;
  return EXIT_SUCCESS;
}

// Get the derivative of the capacitance matrix
int EEQBCModel::get_dcmatdr(
  const TMolecule &mol,  // molecular geometry
  const TIVector &realIdx,  // The real atom indices (for excluding dummy atoms)
  const TMatrix<double> &dist, // atom distances
  TMatrix<double> &dcmatdr  // Out: Capacitance matrix
) {
  int iat;  // atom type of i
  int jat;  // atom type of j
  int ij_min, ij_max; // min and max from i and j
  int ij_at;  // pairwise index for atom types of i and j
  double c_ij;  // Capacitance for bond between atoms i and j
  double dist_ij; // distance between atoms i and j
  double rvdw_ijat;  // pairwise van der Waals radii for ij atom types
  double cap_ij;  // product of bond capacitances for atom types of i and j
  TVector<double> vec; // Vector from i to j
  vec.NewVec(3);
  TVector<double> dcdr_ij; // Part ij of capacitance derivative
  dcdr_ij.NewVec(3);
  const int n_atoms = realIdx.Max() + 1;
  dcmatdr.NewMat(n_atoms, 3 * n_atoms);
  for (int i = 0, ii = 0; i < mol.NAtoms; i++)
  {
    ii = realIdx(i);
    if (ii < 0) continue;
    iat = mol.ATNO(i);
    for (int j = 0, jj = 0; j < i; j++)
    {
      jj = realIdx(j);
      if (jj < 0) continue;
      jat = mol.ATNO(j);
      for (int c = 0; c < 3; c++) {
        vec(c) = mol.CC(jj, c) - mol.CC(ii, c);
      }
      dist_ij = dist(ii,jj);
      ij_min = std::min(iat, jat);
      ij_max = std::max(iat, jat);
      ij_at = ij_min + ij_max * (ij_max - 1)/2 - 1;
      rvdw_ijat = rvdw[ij_at];
      cap_ij = cap[iat] * cap[jat];
      get_dcpair( dist_ij, rvdw_ijat, cap_ij, vec, dcdr_ij);
      for (int c = 0; c < 3; c++) {
      // Calculate Off-diagonal elements; bond capacitances
      dcmatdr(jj, 3*ii+c) = - dcdr_ij(c);
      dcmatdr(ii, 3*jj+c) = + dcdr_ij(c);
      // Calculate diagonal elements; self-capacitance as the negative sum of bond capacitances
      dcmatdr(ii, 3*ii+c) = dcmatdr(ii, 3*ii+c) + dcdr_ij(c);
      dcmatdr(jj, 3*jj+c) = dcmatdr(jj, 3*jj+c) - dcdr_ij(c);
      }
    }
  }

  return EXIT_SUCCESS;
}

// Get the right-hand side vector b
// A * q = b with Coulomb matrix A and partial charges q
// b = X * C with electronegativity vector X and capacitance matrix C
int EEQBCModel::get_xvec(
  const TMolecule &mol,  // molecular geometry
  const TIVector &realIdx,  // The real atom indices (for excluding dummy atoms)
  const TMatrix<double> &dist,  // atom distances
  const TVector<double> &cn,
  TMatrix<double> &cmat, // capacitance matrix
  int charge,  // total charge of the system
  TVector<double> &qloc,
  TVector<double> &xvec  // Out: electronegativity vector
) {
  int info{0};
  const int n_atoms = realIdx.Max() + 1;
  TVector<double> x_tmp;  // dummy for xvec, has dimension N+1 including the constraint
  TMatrix<double> dqlocdr;
  x_tmp.NewVector(n_atoms+1);
  int i_atno;  // atomic number of atom i

  // get local charge
  info = get_qloc(mol, realIdx, dist, charge, qloc, dqlocdr, false);
  if (info != EXIT_SUCCESS) return info;

  for (int i = 0, ii = 0; i < mol.NAtoms; i++)
  {
    ii = realIdx(i);
    if (ii < 0) continue;
    i_atno = mol.ATNO(i);
    x_tmp(ii) = - chi[i_atno] + kcnchi[i_atno]*cn(ii) + kqchi[i_atno]*qloc(ii);
  }
  
  xvec.NewVector(n_atoms + 1);
  xvec(n_atoms) = charge;
  for (int i = 0, ii = 0; i < mol.NAtoms; i++)
  {
    ii = realIdx(i);
    if (ii < 0) continue;
    for (int j = 0, jj = 0; j < mol.NAtoms; j++)
    {
      jj = realIdx(j);
      if (jj < 0) continue;
      xvec(ii) = xvec(ii) + cmat(ii, jj) * x_tmp(jj);
    }
  }

  return EXIT_SUCCESS; 
}

// Get the right-hand side vector b and its derivative w.r.t. CN
// A * q = b with Coulomb matrix A and partial charges q
// b = X * C with electronegativity vector X and capacitance matrix C
int EEQBCModel::get_xvec_derivs(
  const TMolecule &mol,  // molecular geometry
  const TIVector &realIdx,  // The real atom indices (for excluding dummy atoms)
  const TMatrix<double> &dist,  // atom distances
  const TVector<double> &cn,  // coordination number (CN)
  const TMatrix<double> &dcndr, // coordination number derivative
  TMatrix<double> &cmat, // capacitance matrix
  int charge,  // total charge of the system
  TVector<double> &xvec,  // Out: electronegativity vector
  TMatrix<double> &dxvecdr,  // derivative of the entire right hand side w.r.t. atom positions
  TVector<double> &qloc,
  TMatrix<double> &dqlocdr,  // derivative of qloc w.r.t. atom positions
  TMatrix<double> &dcmatdr
) {
  int info{0};
  const int n_atoms = realIdx.Max() + 1;
  int i_atno;  // atomic number of atom i
  TVector<double> x_tmp;  // dummy for xvec, has dimension N+1 including the constraint
  TMatrix<double> dxvecdr_tmp;
  TMatrix<double> cmat_tmp;
  x_tmp.NewVector(n_atoms+1);
  dxvecdr_tmp.NewMat(n_atoms, 3*n_atoms);
  cmat_tmp.NewMat(n_atoms, n_atoms);
  dxvecdr.NewMat(n_atoms, 3*n_atoms);

  // calculate derivative of the capacitance
  info = get_dcmatdr(mol, realIdx, dist, dcmatdr);
  if (info != EXIT_SUCCESS) return info;

  // get local charge
  bool lgrad = true; // calculate dqlocdr
  info = get_qloc(mol, realIdx, dist, charge, qloc, dqlocdr, lgrad);
  if (info != EXIT_SUCCESS) return info;
  for (int i=0, ii = 0; i < mol.NAtoms; i++)
  {
    ii = realIdx(i);
    if (ii < 0) continue;
    i_atno = mol.ATNO(i);
    for (int j = 0, jj = 0; j < mol.NAtoms; j++)
    {
      jj = realIdx(j);
      if (jj < 0) continue;
      // setup C-matrix with correct size
      cmat_tmp(ii,jj) = cmat(ii,jj);
      for (int c = 0; c < 3; c++)
      {
        dxvecdr_tmp(ii, 3*jj+c) = kcnchi[i_atno] * dcndr(ii, 3*jj+c) + dxvecdr_tmp(ii, 3*jj+c);
        dxvecdr_tmp(ii, 3*jj+c) = kqchi[i_atno] * dqlocdr(ii, 3*jj+c) + dxvecdr_tmp(ii, 3*jj+c);
      }
    }
  }
  
  BLAS_Add_Mat_x_Mat(dxvecdr, cmat_tmp, dxvecdr_tmp, false, false, 1.0);

  for (int i = 0, ii = 0; i < mol.NAtoms; i++)
  {
    ii = realIdx(i);
    if (ii < 0) continue;
    i_atno = mol.ATNO(i);
    x_tmp(ii) = - chi[i_atno] + kcnchi[i_atno]*cn(ii) + kqchi[i_atno]*qloc(ii);
  }
  
  xvec.NewVector(n_atoms + 1);
  xvec(n_atoms) = charge;
  BLAS_Add_Mat_x_Vec(xvec, cmat, x_tmp, false, 1.0);

  for (int i = 0, ii = 0; i < mol.NAtoms; i++)
  {
    ii = realIdx(i);
    if (ii < 0) continue;
    for (int j = 0, jj = 0; j < i; j++)
    {
      jj = realIdx(j);
      if (jj < 0) continue;
      for (int c = 0; c < 3; c++) {
        // setup dCij/dR * Xj
        dxvecdr(ii, 3*ii+c) += x_tmp(jj) * dcmatdr(jj, 3*ii+c);
        dxvecdr(jj, 3*jj+c) += x_tmp(ii) * dcmatdr(ii, 3*jj+c);
        dxvecdr(jj, 3*ii+c) += (x_tmp(ii)-x_tmp(jj)) * dcmatdr(jj, 3*ii+c);
        dxvecdr(ii, 3*jj+c) += (x_tmp(jj)-x_tmp(ii)) * dcmatdr(ii, 3*jj+c);
      }
    }
    dxvecdr(ii, 3*ii  ) += x_tmp(ii) * dcmatdr(ii, 3*ii  );
    dxvecdr(ii, 3*ii+1) += x_tmp(ii) * dcmatdr(ii, 3*ii+1);
    dxvecdr(ii, 3*ii+2) += x_tmp(ii) * dcmatdr(ii, 3*ii+2);
  }

  return EXIT_SUCCESS; 
}

} // namespace multicharge
