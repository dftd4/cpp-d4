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
#include <dftd_geometry.h>
#include <dftd_matrix.h>
#include <dftd_model.h>
#include <dftd_ncoord.h>
#include <dftd_readxyz.h>

#include "molecules.h"
#include "test_ghost.h"
#include "util.h"

using namespace dftd4;
using namespace multicharge;

int test_water(
  int n,
  const char atoms[][3],
  const double coord[],
  const double ref_grad[],
  const TIVector &realIdx
) {
  int info;
  multicharge::EEQModel chrg_model;

  // assemble molecule
  int charge{0};
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (info != EXIT_SUCCESS) return info;

  // number of real atoms
  int nat = realIdx.Max() + 1;

  // distances
  TMatrix<double> dist;
  dist.New(nat, nat);
  calc_distances(mol, realIdx, dist);

  // COORDINATION NUMBER CHECK

  // erf-CN without cutoff
  TVector<double> cn;
  TMatrix<double> dcndr;
  NCoordErfD4 ncoord_erf(7.5, 1.0, 9999.9);
  info = ncoord_erf.get_ncoord(mol, realIdx, dist, cn, dcndr, false);
  if (info != EXIT_SUCCESS) return info;

  // compare to ref
  for (int i = 0; i != nat; i++) {
    if (check(cn(i), water_dimer_ref_cn[i]) == EXIT_FAILURE) {
      print_fail("GHOST: CN_D4", cn(i), water_dimer_ref_cn[i]);
      return EXIT_FAILURE;
    }
  }

  ///////////////////////
  // EEQ CHARGES CHECK //
  ///////////////////////

  TVector<double> q;    // partial charges from EEQ model
  TMatrix<double> dqdr; // derivative of partial charges
  q.NewVector(nat);

  // calculate partial charges from EEQ model
  info = chrg_model.get_charges(
    mol, realIdx, dist, charge, cn_eeq_default, q, dqdr, false
  );
  if (info != EXIT_SUCCESS) return info;

  // compare to ref
  for (int i = 0; i != nat; i++) {
    if (check(q(i), water_dimer_ref_q[i]) == EXIT_FAILURE) {
      print_fail("GHOST: EEQ", q(i), water_dimer_ref_q[i]);
      return EXIT_FAILURE;
    }
  }

  ////////////////////////////////
  // TWO-BODY DISPERSION ENERGY //
  ////////////////////////////////

  // two-body dispersion energy
  double ref{-2.3184927839693876E-004};

  double energy{0.0};
  TCutoff cutoff;
  TD4Model d4;

  // BP86 parameters
  dparam par;
  d4par("bp86", par, true);
  par.s9 = 0.0;

  // dispersion main function
  info = get_dispersion(mol, realIdx, charge, d4, par, cutoff, energy, nullptr);
  if (info != EXIT_SUCCESS) return info;

  if (check(energy, ref, 1e-8) == EXIT_FAILURE) {
    print_fail("GHOST: Two-body Energy", energy, ref);
    return EXIT_FAILURE;
  }

  ////////////////////////////
  // FULL DISPERSION ENERGY //
  ////////////////////////////

  energy = 0.0;
  ref = -2.3184926184127263E-004;

  // dispersion main function
  par.s9 = 1.0;
  info = get_dispersion(mol, realIdx, charge, d4, par, cutoff, energy, nullptr);
  if (info != EXIT_SUCCESS) return info;

  if (check(energy, ref, 1e-8) == EXIT_FAILURE) {
    print_fail("GHOST: Full Energy", energy, ref);
    return EXIT_FAILURE;
  }

  /////////////////////////
  // DISPERSION GRADIENT //
  /////////////////////////

  bool lgrad{true};

  // analytical gradient
  double *d4grad = new double[3 * mol.NAtoms];
  for (int i = 0; i < 3 * mol.NAtoms; i++) {
    d4grad[i] = 0.0;
  }

  info = ncoord_erf.get_ncoord(mol, realIdx, dist, cn, dcndr, lgrad);
  if (info != EXIT_SUCCESS) return info;

  dqdr.NewMatrix(3 * nat, nat);
  info = chrg_model.get_charges(
    mol, realIdx, dist, charge, cutoff.cn_eeq, q, dqdr, lgrad
  );
  if (info != EXIT_SUCCESS) return info;

  info = get_dispersion(mol, realIdx, charge, d4, par, cutoff, energy, d4grad);
  if (info != EXIT_SUCCESS) return info;

  for (int i = 0; i < 3 * mol.NAtoms; i++) {
    if (check(d4grad[i], ref_grad[i], 1e-8) == EXIT_FAILURE) {
      print_fail("GHOST: Gradient", d4grad[i], ref_grad[i]);
      delete[] d4grad;
      return EXIT_FAILURE;
    }
  }

  delete[] d4grad;
  return EXIT_SUCCESS;
}

int test_ghost() {
  int info;

  // dummy for masking of ghost atoms
  TIVector realidx;
  realidx.NewVec(water_dimer_n);
  realidx(0) = 0;
  realidx(1) = 1;
  realidx(2) = 2;
  realidx(3) = -1;
  realidx(4) = -1;
  realidx(5) = -1;

  info = test_water(
    water_dimer_n,
    water_dimer_atoms,
    water_dimer_coord,
    water_dimer_ref_grad,
    realidx
  );
  if (info != EXIT_SUCCESS) return info;

  realidx(0) = -1;
  realidx(1) = 0;
  realidx(2) = -1;
  realidx(3) = 1;
  realidx(4) = 2;
  realidx(5) = -1;

  info = test_water(
    water_ghost_n,
    water_ghost_atoms,
    water_ghost_coord,
    water_ghost_ref_grad,
    realidx
  );
  if (info != EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
