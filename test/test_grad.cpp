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
#include <dftd_cutoff.h>
#include <dftd_damping.h>
#include <dftd_dispersion.h>
#include <dftd_eeq.h>
#include <dftd_matrix.h>
#include <dftd_model.h>
#include <dftd_ncoord.h>

#include "molecules.h"
#include "test_grad.h"
#include "util.h"

using namespace dftd4;

int test_numgrad(TMolecule &mol, const int charge, const dparam &par) {
  int info{0};
  double energy{0.0}, er{0.0}, el{0.0};
  double step{1.0e-6};
  double thr{1e-10};

  TMatrix<double> numgrad;
  numgrad.New(mol.NAtoms, 3);

  TCutoff cutoff;
  TD4Model d4;

  // masking (nothing excluded)
  TVector<int> realIdx;
  realIdx.NewVec(mol.NAtoms);
  int nat = 0;
  for (int i = 0; i != mol.NAtoms; i++) {
    realIdx(i) = nat;
    nat++;
  }

  // numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      er = 0.0;
      el = 0.0;

      mol.CC(i, c) += step;
      get_dispersion(mol, realIdx, charge, d4, par, cutoff, er, nullptr);

      mol.CC(i, c) = mol.CC(i, c) - 2 * step;
      get_dispersion(mol, realIdx, charge, d4, par, cutoff, el, nullptr);

      mol.CC(i, c) = mol.CC(i, c) + step;
      numgrad(i, c) = 0.5 * (er - el) / step;
    }
  }

  // analytical gradient
  double *d4grad = new double[3 * mol.NAtoms];
  for (int i = 0; i < 3 * mol.NAtoms; i++) {
    d4grad[i] = 0.0;
  }
  info = get_dispersion(mol, realIdx, charge, d4, par, cutoff, energy, d4grad);
  if (info != EXIT_SUCCESS) return info;

  // check translational invariance of analytical gradient
  if (is_trans_invar(mol, d4grad) != EXIT_SUCCESS) return EXIT_FAILURE;

  // compare against numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      if (check(d4grad[3 * i + c], numgrad(i, c), thr) != EXIT_SUCCESS) {
        print_fail("Gradient mismatch", d4grad[3 * i + c], numgrad(i, c));
        delete[] d4grad;
        return EXIT_FAILURE;
      }
    }
  }

  delete[] d4grad;

  return EXIT_SUCCESS;
}

int is_trans_invar(const TMolecule &mol, double gradient[]) {
  double xsum{0.0};
  double ysum{0.0};
  double zsum{0.0};
  for (int i = 0; i < mol.NAtoms; i++) {
    xsum += gradient[3 * i];
    ysum += gradient[3 * i + 1];
    zsum += gradient[3 * i + 2];
  }

  if (check(xsum, 0.0) != EXIT_SUCCESS) {
    print_fail("Translational invariance (x)", xsum, 0.0);
    return EXIT_FAILURE;
  }
  if (check(ysum, 0.0) != EXIT_SUCCESS) {
    print_fail("Translational invariance (y)", ysum, 0.0);
    return EXIT_FAILURE;
  }
  if (check(zsum, 0.0) != EXIT_SUCCESS) {
    print_fail("Translational invariance (z)", zsum, 0.0);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int test_pbed4_mb01() {
  // PBE-D4(EEQ) parameters
  dparam par;
  par.s6 = 1.0;
  par.s9 = 0.0; // no MBD!
  par.alp = 16.0;
  par.s8 = 0.95948085;
  par.a1 = 0.38574991;
  par.a2 = 4.80688534;
  par.s10 = 0.0;

  // assemble molecule
  int charge = mb16_43_01_charge;
  TMolecule mol;
  int info =
    get_molecule(mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord, mol);
  if (info != EXIT_SUCCESS) return info;

  return test_numgrad(mol, charge, par);
};

int test_bp86d4atm_water() {
  // PBE-D4(EEQ)-ATM parameters
  dparam par;
  d4par("bp86", par, true);

  // assemble molecule
  int charge = water_charge;
  TMolecule mol;
  int info = get_molecule(water_n, water_atoms, water_coord, mol);
  if (info != EXIT_SUCCESS) return info;

  return test_numgrad(mol, charge, par);
}

int test_tpss0d4mbd_rost61m1() {
  // TPSS0-D4(EEQ)-MBD parameters
  dparam par;
  d4par("tpss0", par, false);

  // assemble molecule
  int charge = rost61_m1_charge;
  TMolecule mol;
  int info = get_molecule(rost61_m1_n, rost61_m1_atoms, rost61_m1_coord, mol);
  if (info != EXIT_SUCCESS) return info;

  return test_numgrad(mol, charge, par);
}

int test_numgrad_dqdr_eeq(int n, const char atoms[][3], const double coord[]) {
  // assemble molecule
  int info;
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (info != EXIT_SUCCESS) return info;
  TVector<double> q_r, q_l;
  double step{1.0e-6}; // accurate up to 1.0E-8
  double thr{1.2e-8};

  TMatrix<double> dist;
  dist.NewMat(mol.NAtoms, mol.NAtoms);
  TMatrix<double> num_dqdr; // numerical gradient of the partial charges
  num_dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);
  TMatrix<double> analytic_dqdr; // analytical gradient of the partial charges
  analytic_dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);
  multicharge::EEQModel eeq_model;
  TVector<double> q;
  TMatrix<double> dqdr;
  q.NewVec(mol.NAtoms);
  dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);

  TCutoff cutoff;

  // masking (nothing excluded)
  TVector<int> realIdx;
  realIdx.NewVec(mol.NAtoms);
  int nat = 0;
  for (int i = 0; i != mol.NAtoms; i++) {
    realIdx(i) = nat;
    nat++;
  }

  // analytical gradient
  calc_distances(mol, realIdx, dist);
  info = eeq_model.get_charges(
    mol, realIdx, dist, 0, cutoff.cn_eeq, q, analytic_dqdr, true
  );
  if (info != EXIT_SUCCESS) return info;

  // calculate numerical gradient via finite difference method
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      // calculate forward point
      q_r.NewVec(n);
      mol.CC(i, c) += step;
      calc_distances(mol, realIdx, dist);
      eeq_model.get_charges(
        mol, realIdx, dist, 0, cutoff.cn_eeq, q_r, dqdr, false
      );

      // calculate backward point
      q_l.NewVec(n);
      mol.CC(i, c) = mol.CC(i, c) - 2 * step;
      calc_distances(mol, realIdx, dist);
      eeq_model.get_charges(
        mol, realIdx, dist, 0, cutoff.cn_eeq, q_l, dqdr, false
      );

      // calculate numerical gradient as finite difference
      mol.CC(i, c) = mol.CC(i, c) + step;
      for (int j = 0; j < mol.NAtoms; j++) {
        num_dqdr(3 * i + c, j) = 0.5 * (q_r(j) - q_l(j)) / step;
      }
    }
  }

  // compare against numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      for (int j = 0; j < mol.NAtoms; j++) {
        if (check(analytic_dqdr(3 * i + c, j), num_dqdr(3 * i + c, j), thr) !=
            EXIT_SUCCESS) {
          print_fail(
            "Gradient mismatch for dqdr in EEQ.\n",
            analytic_dqdr(3 * i + c, j),
            num_dqdr(3 * i + c, j)
          );
          return EXIT_FAILURE;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

int test_numgrad_dqdr_eeqbc(
  int n,
  const char atoms[][3],
  const double coord[]
) {
  // assemble molecule
  int info;
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (info != EXIT_SUCCESS) return info;
  TVector<double> q_r, q_l;
  double step{1.0e-6}; // accurate up to 1.0E-8
  double thr{4.0e-8};

  TMatrix<double> dist;
  dist.NewMat(mol.NAtoms, mol.NAtoms);
  TMatrix<double> num_dqdr; // numerical gradient of the partial charges
  num_dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);
  TMatrix<double> analytic_dqdr; // analytical gradient of the partial charges
  analytic_dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);
  multicharge::EEQBCModel eeqbc_model;
  TVector<double> q;
  TMatrix<double> dqdr;
  q.NewVec(mol.NAtoms);
  dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);

  TCutoff cutoff;

  // masking (nothing excluded)
  TVector<int> realIdx;
  realIdx.NewVec(mol.NAtoms);
  int nat = 0;
  for (int i = 0; i != mol.NAtoms; i++) {
    realIdx(i) = nat;
    nat++;
  }

  // analytical gradient
  calc_distances(mol, realIdx, dist);
  info = eeqbc_model.get_charges(
    mol, realIdx, dist, 0, cutoff.cn_eeq, q, analytic_dqdr, true
  );
  if (info != EXIT_SUCCESS) return info;

  // calculate numerical gradient via finite difference method
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      // calculate forward point
      q_r.NewVec(n);
      mol.CC(i, c) += step;
      calc_distances(mol, realIdx, dist);
      eeqbc_model.get_charges(
        mol, realIdx, dist, 0, cutoff.cn_eeq, q_r, dqdr, false
      );

      // calculate backward point
      q_l.NewVec(n);
      mol.CC(i, c) = mol.CC(i, c) - 2 * step;
      calc_distances(mol, realIdx, dist);
      eeqbc_model.get_charges(
        mol, realIdx, dist, 0, cutoff.cn_eeq, q_l, dqdr, false
      );

      // calculate numerical gradient as finite difference
      mol.CC(i, c) = mol.CC(i, c) + step;
      for (int j = 0; j < mol.NAtoms; j++) {
        num_dqdr(3 * i + c, j) = 0.5 * (q_r(j) - q_l(j)) / step;
      }
    }
  }

  // compare against numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      for (int j = 0; j < mol.NAtoms; j++) {
        if (check(analytic_dqdr(3 * i + c, j), num_dqdr(3 * i + c, j), thr) !=
            EXIT_SUCCESS) {
          print_fail(
            "Gradient mismatch for dqdr in EEQ-BC.\n",
            analytic_dqdr(3 * i + c, j),
            num_dqdr(3 * i + c, j)
          );
          return EXIT_FAILURE;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

int test_grad() {
  int info{0};

  info =
    test_numgrad_dqdr_eeq(mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord);
  if (info != EXIT_SUCCESS) return info;

  info =
    test_numgrad_dqdr_eeqbc(mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord);
  if (info != EXIT_SUCCESS) return info;

  info = test_pbed4_mb01();
  if (info != EXIT_SUCCESS) return info;

  info = test_bp86d4atm_water();
  if (info != EXIT_SUCCESS) return info;

  info = test_tpss0d4mbd_rost61m1();
  if (info != EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
};
