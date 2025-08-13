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

  // get correct cn for current charge model
  get_cn(mol, realIdx, dist, lgrad);

  // corresponds to model%solve in Fortran
  info =
    eeq_chrgeq(mol, realIdx, dist, charge, q, dqdr, lgrad, lverbose);
  if (info != EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
};

int ChargeModel::eeq_chrgeq(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const int &charge,
  TVector<double> &qvec,
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

  info = get_vrhs(mol, realIdx, charge, dist, xvec, dxdcn, lgrad);
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

        dAmat(3 * jj,     ii) -= ncoord_erf.dcndr(ii, 3 * jj    ) * dxdcn(ii);
        dAmat(3 * jj + 1, ii) -= ncoord_erf.dcndr(ii, 3 * jj + 1) * dxdcn(ii);
        dAmat(3 * jj + 2, ii) -= ncoord_erf.dcndr(ii, 3 * jj + 2) * dxdcn(ii);
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
  TVector<double> &xvec,
  TVector<double> &dxvec,
  bool lgrad
) {
  double tmp{0.0};
  int izp;
  int nat = realIdx.Max() + 1;

  if (lgrad) {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(ncoord_erf.cn(ii) + small);
      xvec(ii) = -xi[izp] + tmp * ncoord_erf.cn(ii);
      dxvec(ii) = 0.5 * tmp;
    }
    dxvec(nat) = 0.0;
  } else {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(ncoord_erf.cn(ii) + small);
      xvec(ii) = -xi[izp] + tmp * ncoord_erf.cn(ii);
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

      dAmat(3 * ii,     jj) =  dgx * q(ii);
      dAmat(3 * ii + 1, jj) =  dgy * q(ii);
      dAmat(3 * ii + 2, jj) =  dgz * q(ii);
      dAmat(3 * jj,     ii) = -dgx * q(jj);
      dAmat(3 * jj + 1, ii) = -dgy * q(jj);
      dAmat(3 * jj + 2, ii) = -dgz * q(jj);
    }
  }

  return EXIT_SUCCESS;
};

// Calculate the coordination number, forwarding to get_ncoord
NCoordErf& EEQModel::get_cn(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  bool lgrad

) 
{
  int info;

  // get the EEQ coordination number
  info = ncoord_erf.get_ncoord(mol, realIdx, dist, lgrad);
  if (info == EXIT_SUCCESS) {
    return ncoord_erf;
  } else {
    throw std::runtime_error("EEQModel::get_cn: Failed to compute coordination numbers.");
  }
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
  ncoord_erf = NCoordErf(2.0, 0.75, 25.0, 1.0, -1.0, multicharge_param::eeqbc::eeqbc_cov_radii);
}

// Calculate the right-hand side (RHS) vector of the system of linear equations for EEQ-BC
int EEQBCModel::get_vrhs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int &charge,
  const TMatrix<double> &dist,
  TVector<double> &xvec,
  TVector<double> &dxvec,
  bool lgrad
) {
  int info;
  int n_atoms = realIdx.Max() + 1;
  // calculate the capacitance matrix
  cmat.NewMatrix(n_atoms + 1, n_atoms + 1);
  info = get_cmat(mol, realIdx, dist, cmat);
  if (info != EXIT_SUCCESS) {
     printf("EEQBCModel::get_vrhs: Failed to calculate the capacitance matrix.");
     return info;
  }
  
  // calculate the right-hand side
  info = get_xvec(mol, realIdx, dist, cmat, ncoord_erf, charge, xvec);
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
  TMatrix<double> &Amat
) const {
  int nat = mol.NAtoms;
  int iat, jat; // atomic numbers
  double norm_cn; // coordination number normalization factor
  double r, radi, radj, gamij2, tmp;


  for (int i = 0; i < nat; i++) {
    //
    iat = mol.ATNO(i);
    norm_cn = 1.0 / pow(avg_cn[iat], norm_exp);
    radi = rad[iat] * (1.0 - kcnrad*ncoord_erf.cn(i)*norm_cn);

    for (int j = 0; j < i; j++) {
      //
      jat = mol.ATNO(j);
      r = dist(i, j);
      norm_cn =  1.0 / pow(avg_cn[jat], norm_exp);
      radj = rad[jat] * (1.0 - kcnrad*ncoord_erf.cn(j)*norm_cn);
      gamij2 = 1.0 / (radi * radi + radj * radj);
      tmp = cmat(j, i) * erf(r * sqrt(gamij2)) / r;
      Amat(j, i) = Amat(j, i) + tmp;
      Amat(i, j) = Amat(i, j) + tmp;
    }
    tmp = eta[iat] + kqeta[iat] * qloc(i) + sqrt(2.0/pi) / radi;
    Amat(i, i) = Amat(i, i) + tmp * cmat(i, i) + 1.0;
    // Set entries for charge constraint
    Amat(nat, i) = 1.0;
    Amat(i, nat) = 1.0;
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
  TMatrix<double> &dAmat,
  TMatrix<double> &atrace
) const {
  
  return EXIT_SUCCESS;
};

// Get purely geometry-dependent local charges
int EEQBCModel::get_qloc(
  const TMolecule &mol,
  const TMatrix<double> &dist,
  const double q_tot,  // total system charge
  TVector<double> &qloc
) 
{
  const double cutoff = 25.0;
  bool lgrad = false;
  NCoordErfEN ncoord_erf_en(2.0, 0.75, 25.0, -1.0, -1.0);
  qloc.NewVector(mol.NAtoms);
  const double q_tot_norm = q_tot/mol.NAtoms;

  ncoord_erf_en.get_ncoord(mol, dist, cutoff, lgrad);

  for (int i=0; i < mol.NAtoms; i++) {
    qloc(i) = ncoord_erf_en.cn(i) + q_tot_norm;
  }

  return EXIT_SUCCESS;
};


// Calculate the coordination number, forwarding to get_ncoord
NCoordErf& EEQBCModel::get_cn(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  bool lgrad
) 
{
  int info;

  info = ncoord_erf.get_ncoord(mol, realIdx, dist, lgrad);
  if (info == EXIT_SUCCESS) {
    return ncoord_erf;
  } else {
    throw std::runtime_error("EEQModel::get_cn: Failed to compute coordination numbers.");
  }
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

// Get the capacitance matrix
int EEQBCModel::get_cmat(
  const TMolecule &mol,  // molecular geometry
  const TIVector &realIdx,
  const TMatrix<double> &dist, // atom distances
  TMatrix<double> &cmat  // Out: Capacitance matrix
) {
  int iat;  // atom type of i
  int jat;  // atom type of j
  double c_ij;  // Capacitance for bond between atoms i and j
  double dist_ij; // distance between atoms i and j
  const int n_atoms = realIdx.Max() + 1;
  for (int i = 0; i < n_atoms; i++) 
  {
    iat = mol.ATNO(i);
    for (int j = 0; j < i; j++)
    {
      jat = mol.ATNO(j);
      dist_ij = dist(i,j);
      get_cpair(iat, jat, dist_ij, c_ij);
      // Calulate Off-diagonal elements; bond capacitances
      cmat(i, j) = - c_ij;
      cmat(j, i) = - c_ij;
      // Calculate diagonal elements; self-capacitance as the negative sum of bond capacitances
      cmat(i, i) = cmat(i, i) + c_ij;
      cmat(j, j) = cmat(j, j) + c_ij;
    }
  }
  cmat(n_atoms, n_atoms) = 1.0;
  return EXIT_SUCCESS;
}

// Get the right-hand side vector b
// A * q = b with Coulomb matrix A and partial charges q
// b = X * C with electronegativity vector X and capacitance matrix C
int EEQBCModel::get_xvec(
  const TMolecule &mol,  // molecular geometry
  const TIVector &realIdx,
  const TMatrix<double> &dist,  // atom distances
  TMatrix<double> &cmat, // capacitance matrix
  NCoordErf &ncoord_erf,  // erf()-based coordination number object for EEQ-BC
  int charge,  // total charge of the system
  TVector<double> &xvec  // Out: electronegativity vector
) {
  int info;
  const int n_atoms = realIdx.Max() + 1;
  TVector<double> x_tmp;  // dummy for xvec, has dimension N+1 including the constraint
  x_tmp.NewVector(n_atoms+1);
  int i_atno;  // atomic number of atom i

  // get local charge
  info = get_qloc(mol, dist, charge, qloc);

  for (int i=0; i < n_atoms; i++)
  {
    i_atno = mol.ATNO(i);
    x_tmp(i) = - chi[i_atno] + kcnchi[i_atno]*ncoord_erf.cn(i) + kqchi[i_atno]*qloc(i);
  }
  
  xvec.NewVector(n_atoms + 1);
  xvec(n_atoms) = charge;
  for (int i = 0; i < mol.NAtoms; i++)
  {
    for (int j = 0; j < mol.NAtoms; j++)
    {
      xvec(i) = xvec(i) + cmat(i, j) * x_tmp(j);
    }
  }

  return EXIT_SUCCESS; 
}

// numerical gradient of partial charges w.r.t. atom positions
int EEQBCModel::num_grad_dqdr(
  TMolecule &mol,  // molecular geometry
  const TIVector &realIdx,
  int charge, // total charge of the system
  TMatrix<double> &num_dqdr // numerical gradient
) {
  TVector<double> q_r, q_l;
  double step{1.0e-6};
  int nat = realIdx.Max() +1;
  TMatrix<double> dqdr;
  num_dqdr.NewMat(3 * nat, nat);
  TMatrix<double> dist;
  TVector<double> q_dum;
  for (int i = 0; i < nat; i++) {
    for (int c = 0; c < 3; c++) {
      multicharge::EEQBCModel eeqbc_model_l;
      multicharge::EEQBCModel eeqbc_model_r;
      q_r.NewVec(nat);
      q_l.NewVec(nat);
      q_dum.NewVec(nat);
      dqdr.NewMat(3 * nat, nat);
      mol.CC(i, c) += step;
      dist.NewMatrix(nat, nat);
      calc_distances(mol, realIdx, dist);
      eeqbc_model_r.get_cn(mol, realIdx, dist, false);
      eeqbc_model_r.eeq_chrgeq(mol, realIdx, dist, charge, q_dum, dqdr, false, false);

      for (int k = 0; k < nat; k++) {
        q_r(k) = q_dum(k);
      }

      mol.CC(i, c) = mol.CC(i, c) - 2 * step;
      dist.NewMatrix(nat, nat);
      calc_distances(mol, realIdx, dist);
      q_dum.NewVec(nat);
      dqdr.NewMat(3 * nat, nat);
      eeqbc_model_l.get_cn(mol, realIdx, dist, false);
      eeqbc_model_l.eeq_chrgeq(mol, realIdx, dist, charge, q_dum, dqdr, false, false);
      for (int k = 0; k < nat; k++) {
        q_l(k) = q_dum(k);
      }

      mol.CC(i, c) = mol.CC(i, c) + step;
      for (int j = 0; j < nat; j++) {
        num_dqdr(3 * i + c, j) = 0.5 * (q_r(j) - q_l(j)) / step;

      }
    }
  }
  return EXIT_SUCCESS;
}

} // namespace multicharge
