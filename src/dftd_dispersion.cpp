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

#include "damping/dftd_atm.h"
#include "damping/dftd_rational.h"
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
  const TD4Model &d4,
  const dparam &par,
  const TCutoff cutoff,
  double &energy,
  double *GRAD
) {
  TVector<int> realIdx;
  initializeRealIdx(mol.NAtoms, realIdx);

  return get_dispersion(mol, realIdx, charge, d4, par, cutoff, energy, GRAD);
}

int get_dispersion(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int charge,
  const TD4Model &d4,
  const dparam &par,
  const TCutoff cutoff,
  double &energy,
  double *GRAD
) {

  int info{0};

  int nat = realIdx.Max() + 1;

  TRVector energies; // atom-wise energies
  energies.NewVector(nat);

  info = get_dispersion(mol, realIdx, charge, d4, par, cutoff, energies, GRAD);
  if (info != EXIT_SUCCESS) return info;

  // sum up atom-wise energies
  for (int i = 0; i != nat; i++) {
    energy += energies(i);
  }

  energies.DelVec();

  return EXIT_SUCCESS;
}

int get_dispersion(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int charge,
  const TD4Model &d4,
  const dparam &par,
  const TCutoff cutoff,
  TRVector &energies,
  double *GRAD
) {
  // setup variables
  int info{0};
  bool lmbd = (par.s9 != 0.0);
  bool lgrad = !!GRAD;

  int nat = realIdx.Max() + 1;

  // distances
  TMatrix<double> dist;
  dist.NewMatrix(nat, nat);
  calc_distances(mol, realIdx, dist);

  TVector<double> q;        // partial charges from EEQ model
  TMatrix<double> dqdr;     // derivative of partial charges
  TVector<double> gradient; // derivative of dispersion energy
  NCoordErfD4
    ncoord_erf_d4;       // instance of erf() based coordination number for D4
  TVector<double> cn;    // coordination number
  TMatrix<double> dcndr; // derivative of the coordination number
  multicharge::EEQModel chrg_model; // EEQ charge model

  q.NewVector(nat);
  if (lgrad) {
    dqdr.NewMatrix(3 * nat, nat);
    gradient.NewVector(3 * nat);
  }

  // calculate partial charges from EEQ model
  info = chrg_model.get_charges(
    mol, realIdx, dist, charge, cutoff.cn_eeq, q, dqdr, lgrad
  );
  if (info != EXIT_SUCCESS) return info;

  // get the D4 coordination number
  info = ncoord_erf_d4.get_ncoord(mol, realIdx, dist, cn, dcndr, lgrad);
  if (info != EXIT_SUCCESS) return info;

  // maximum number of reference systems
  int mref{0};
  info = get_max_ref(mol, mref);
  if (info != EXIT_SUCCESS) return info;

  // reference charges
  TMatrix<double> refq;
  refq.NewMat(mref, nat);
  info = d4.set_refq_eeq(mol, realIdx, refq);
  if (info != EXIT_SUCCESS) return info;

  TMatrix<double> gwvec;
  TMatrix<double> dgwdcn;
  TMatrix<double> dgwdq;
  gwvec.NewMatrix(mref, nat);
  if (lgrad) {
    dgwdcn.NewMatrix(mref, nat);
    dgwdq.NewMatrix(mref, nat);
  }
  info = d4.weight_references(
    mol, realIdx, cn, q, refq, gwvec, dgwdcn, dgwdq, lgrad
  );
  if (info != EXIT_SUCCESS) return info;

  TMatrix<double> c6;
  TMatrix<double> dc6dcn;
  TMatrix<double> dc6dq;
  c6.NewMatrix(nat, nat);
  if (lgrad) {
    dc6dcn.NewMatrix(nat, nat);
    dc6dq.NewMatrix(nat, nat);
  }

  info = d4.get_atomic_c6(
    mol, realIdx, gwvec, dgwdcn, dgwdq, c6, dc6dcn, dc6dq, lgrad
  );
  if (info != EXIT_SUCCESS) return info;

  // --------------------------
  // Two-body dispersion energy
  // --------------------------

  TVector<double> dEdcn;
  TVector<double> dEdq;
  if (lgrad) {
    dEdcn.NewVector(nat);
    dEdq.NewVector(nat);
  }

  info = get_dispersion2(
    mol,
    realIdx,
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
  if (info != EXIT_SUCCESS) return info;

  // blas function have no return in orca
  if (lgrad) BLAS_Add_Mat_x_Vec(gradient, dqdr, dEdq, false, 1.0);
  dqdr.DelMat();

  // ----------------------------
  // Three-body dispersion energy
  // ----------------------------

  if (lmbd) {
    // Three-body term is independent of charges
    for (int i = 0; i != nat; i++) {
      q(i) = 0.0;
    }

    // calculate weight references
    gwvec.NewMatrix(mref, nat);
    if (lgrad) {
      dgwdcn.NewMatrix(mref, nat);
      dgwdq.NewMatrix(mref, nat);
    }
    info = d4.weight_references(
      mol, realIdx, cn, q, refq, gwvec, dgwdcn, dgwdq, lgrad
    );
    if (info != EXIT_SUCCESS) return info;

    q.Delete();
    refq.Delete();

    // calculate reference C6 coefficients
    c6.NewMatrix(nat, nat);
    if (lgrad) {
      dc6dcn.NewMatrix(nat, nat);
      dc6dq.NewMatrix(nat, nat);
    }
    info = d4.get_atomic_c6(
      mol, realIdx, gwvec, dgwdcn, dgwdq, c6, dc6dcn, dc6dq, lgrad
    );
    if (info != EXIT_SUCCESS) return info;

    gwvec.DelMat();
    dgwdcn.DelMat();
    dgwdq.DelMat();

    // calculate three-body dispersion
    info = get_dispersion3(
      mol,
      realIdx,
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
    if (info != EXIT_SUCCESS) return info;
  } else {
    q.Delete();
    refq.Delete();
    gwvec.Delete();
    dgwdcn.Delete();
    dgwdq.Delete();
  }

  dist.DelMat();
  c6.DelMat();
  dc6dcn.DelMat();
  dc6dq.DelMat();

  if (lgrad) BLAS_Add_Mat_x_Vec(gradient, dcndr, dEdcn, true, 1.0);

  dEdcn.DelVec();
  dEdq.DelVec();

  // write to input gradient
  if (lgrad) {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      GRAD[3 * i] = gradient(3 * ii);
      GRAD[3 * i + 1] = gradient(3 * ii + 1);
      GRAD[3 * i + 2] = gradient(3 * ii + 2);
    }

    gradient.DelVec();
  }

  return EXIT_SUCCESS;
}

} // namespace dftd4
