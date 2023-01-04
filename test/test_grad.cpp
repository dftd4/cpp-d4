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
#include <dftd_matrix.h>

#include "molecules.h"
#include "test_grad.h"
#include "util.h"

using namespace dftd4;

int test_numgrad(const TMolInfo &dat, TMolecule &mol, const dparam &par) {
  int info{0};
  double energy{0.0}, er{0.0}, el{0.0};
  double step{1.0e-6};
  double thr{1e-10};

  TMatrix<double> numgrad;
  numgrad.New(mol.NAtoms, 3);

  TCutoff cutoff;

  // numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      er = 0.0;
      el = 0.0;

      mol.xyz(i, c) += step;
      get_dispersion(dat, mol, par, cutoff, er, nullptr);

      mol.xyz(i, c) = mol.xyz(i, c) - 2*step;
      get_dispersion(dat, mol, par, cutoff, el, nullptr);

      mol.xyz(i, c) = mol.xyz(i, c) + step;
      numgrad(i, c) = 0.5 * (er - el) / step;
    }    
  }

  // analytical gradient
  double* d4grad = new double[3*mol.NAtoms];
  for (int i = 0; i < 3*mol.NAtoms; i++) d4grad[i] = 0.0;
  info = get_dispersion(dat, mol, par, cutoff, energy, d4grad);
  if (!info == EXIT_SUCCESS) return info;

  // check translational invariance of analytical gradient
  if (is_trans_invar(mol, d4grad) != EXIT_SUCCESS) return EXIT_FAILURE;

  // compare against numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      if (check(d4grad[3*i + c], numgrad(i, c), thr) != EXIT_SUCCESS) {
        print_fail("Gradient mismatch", d4grad[3*i + c], numgrad(i, c));
        return EXIT_FAILURE;
      }
    }
  }

  return EXIT_SUCCESS;
}

int is_trans_invar(const TMolecule& mol, double gradient[]) {
  double xsum{0.0};
  double ysum{0.0};
  double zsum{0.0};
  for (int i = 0; i < mol.NAtoms; i++) {
    xsum += gradient[3*i];
    ysum += gradient[3*i + 1];
    zsum += gradient[3*i + 2];
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

  // assemble molecule
  TMolInfo dat = TMolInfo(mb16_43_01_charge);
  TMolecule mol;
  int info = get_molecule(mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord, mol);
  if (!info == EXIT_SUCCESS) return info;

  return test_numgrad(dat, mol, par);
};

int test_bp86d4atm_water() {
  // PBE-D4(EEQ)-ATM parameters
  dparam par;
  d4par("bp86", par, true);

  // assemble molecule
  TMolInfo dat = TMolInfo(water_charge);
  TMolecule mol;
  int info = get_molecule(water_n, water_atoms, water_coord, mol);
  if (!info == EXIT_SUCCESS) return info;

  return test_numgrad(dat, mol, par);
}

int test_tpss0d4mbd_rost61m1() {
  // TPSS0-D4(EEQ)-MBD parameters
  dparam par;
  d4par("tpss0", par, false);

  // assemble molecule
  TMolInfo dat = TMolInfo(rost61_m1_charge);
  TMolecule mol;
  int info = get_molecule(rost61_m1_n, rost61_m1_atoms, rost61_m1_coord, mol);
  if (!info == EXIT_SUCCESS) return info;

  return test_numgrad(dat, mol, par);
}

int test_grad() {
  int info{0};

  info = test_pbed4_mb01();
  if (!info == EXIT_SUCCESS) return info;

  info = test_bp86d4atm_water();
  if (!info == EXIT_SUCCESS) return info;

  info = test_tpss0d4mbd_rost61m1();
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
};
