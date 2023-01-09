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
 * D4(EEQ)-ATM implementation
 */
#include <cmath>

#include "damping/atm.h"
#include "damping/rational.h"
#include "dftd_cblas.h"
#include "dftd_dispersion.h"
#include "dftd_eeq.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_model.h"
#include "dftd_ncoord.h"
#include "dftd_parameters.h"

namespace dftd4 {

int get_dispersion(
  const TMolecule &mol,
  const int charge,
  const dparam &par,
  TCutoff cutoff,
  double &energy,
  double *GRAD
) {
  // setup variables
  bool lmbd = (par.s9 != 0.0);
  bool lgrad = !!GRAD;
  int info = 0;

  // distances
  TMatrix<double> dist;
  dist.NewMat(mol.NAtoms, mol.NAtoms);
  calc_distances(mol, dist);

  TVector<double> cn;       // D4 coordination number
  TVector<double> q;        // partial charges from EEQ model
  TMatrix<double> dcndr;    // derivative of D4-CN
  TMatrix<double> dqdr;     // derivative of partial charges
  TVector<double> gradient; // derivative of dispersion energy

  cn.NewVec(mol.NAtoms);
  q.NewVec(mol.NAtoms);
  if (lgrad) {
    dcndr.NewMat(3 * mol.NAtoms, mol.NAtoms);
    dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);
    gradient.NewVec(3 * mol.NAtoms);
  }

  // calculate partial charges from EEQ model
  info = get_charges(mol, dist, charge, cutoff.cn_eeq, q, dqdr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // get the D4 coordination number
  info = get_ncoord_d4(mol, dist, cutoff.cn, cn, dcndr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // maximum number of reference systems
  int mref{0};
  info = get_max_ref(mol, mref);
  if (!info == EXIT_SUCCESS) return info;

  TMatrix<double> gwvec;
  TMatrix<double> dgwdcn;
  TMatrix<double> dgwdq;
  gwvec.NewMat(mref, mol.NAtoms);
  if (lgrad) {
    dgwdcn.NewMat(mref, mol.NAtoms);
    dgwdq.NewMat(mref, mol.NAtoms);
  }
  info = weight_references(mol, cn, q, gwvec, dgwdcn, dgwdq, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  TMatrix<double> c6;
  TMatrix<double> dc6dcn;
  TMatrix<double> dc6dq;
  c6.NewMat(mol.NAtoms, mol.NAtoms);
  if (lgrad) {
    dc6dcn.NewMat(mol.NAtoms, mol.NAtoms);
    dc6dq.NewMat(mol.NAtoms, mol.NAtoms);
  }

  info = get_atomic_c6(mol, gwvec, dgwdcn, dgwdq, c6, dc6dcn, dc6dq, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // --------------------------
  // Two-body dispersion energy
  // --------------------------

  TVector<double> dEdcn;
  TVector<double> dEdq;
  TVector<double> energies;
  energies.NewVec(mol.NAtoms);
  if (lgrad) {
    dEdcn.NewVec(mol.NAtoms);
    dEdq.NewVec(mol.NAtoms);
  }

  info = get_dispersion2(
    mol,
    dist,
    cutoff.disp2,
    par,
    c6,
    dc6dcn,
    dc6dq,
    energies,
    dEdcn,
    dEdq,
    gradient,
    lgrad
  );
  if (!info == EXIT_SUCCESS) return info;

  if (lgrad) {
    info = BLAS_Add_Mat_x_Vec(gradient, dqdr, dEdq, false, 1.0);
    if (!info == EXIT_SUCCESS) return info;
  }

  dqdr.Delete();

  // ----------------------------
  // Three-body dispersion energy
  // ----------------------------

  if (lmbd) {
    // Three-body term is independent of charges
    for (int i = 0; i != mol.NAtoms; i++) {
      q(i) = 0.0;
    }

    // calculate weight references
    gwvec.NewMat(mref, mol.NAtoms);
    if (lgrad) {
      dgwdcn.NewMat(mref, mol.NAtoms);
      dgwdq.NewMat(mref, mol.NAtoms);
    }
    info = weight_references(mol, cn, q, gwvec, dgwdcn, dgwdq, lgrad);
    if (!info == EXIT_SUCCESS) return info;

    cn.Delete();
    q.Delete();

    // calculate reference C6 coefficients
    c6.NewMat(mol.NAtoms, mol.NAtoms);
    if (lgrad) {
      dc6dcn.NewMat(mol.NAtoms, mol.NAtoms);
      dc6dq.NewMat(mol.NAtoms, mol.NAtoms);
    }
    info = get_atomic_c6(mol, gwvec, dgwdcn, dgwdq, c6, dc6dcn, dc6dq, lgrad);
    if (!info == EXIT_SUCCESS) return info;

    gwvec.Delete();
    dgwdcn.Delete();
    dgwdq.Delete();

    // calculate three-body dispersion
    info = get_dispersion3(
      mol,
      dist,
      cutoff.disp3,
      par,
      c6,
      dc6dcn,
      dc6dq,
      energies,
      dEdcn,
      dEdq,
      gradient,
      lgrad
    );
    if (!info == EXIT_SUCCESS) return info;
  } else {
    cn.Delete();
    q.Delete();
    gwvec.Delete();
    dgwdcn.Delete();
    dgwdq.Delete();
  }

  dist.Delete();
  c6.Delete();
  dc6dcn.Delete();
  dc6dq.Delete();

  if (lgrad) {
    info = BLAS_Add_Mat_x_Vec(gradient, dcndr, dEdcn, false, 1.0);
    if (!info == EXIT_SUCCESS) return info;
  }

  dcndr.Delete();
  dEdcn.Delete();
  dEdq.Delete();

  // sum up atom-wise energies
  for (int i = 0; i != mol.NAtoms; i++) {
    energy += energies(i);
  }
  energies.Delete();

  // write to input gradient
  if (lgrad) {
    for (int i = 0; i != 3 * mol.NAtoms; i++) {
      GRAD[i] = gradient(i);
    }
    gradient.Delete();
  }

  return EXIT_SUCCESS;
}

} // namespace dftd4
