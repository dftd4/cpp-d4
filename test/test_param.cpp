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
#include <cmath>
#include <limits>
#include <string>

#include <dftd_cutoff.h>
#include <dftd_damping.h>
#include <dftd_geometry.h>
#include <dftd_model.h>
#include <dftd_parameters.h>

#include "molecules.h"
#include "test_param.h"
#include "util.h"

using namespace dftd4;

int test_rational_damping(const double ref[], TCutoff cutoff) {
  int info;
  double energy;
  dparam par;

  // assemble molecule
  int charge = upu23_0a_charge;
  TMolecule mol;
  info = get_molecule(upu23_0a_n, upu23_0a_atoms, upu23_0a_coord, mol);
  if (info != EXIT_SUCCESS) return info;

  // create default D4 model
  TD4Model d4;

  // masking (nothing excluded)
  TVector<int> realIdx;
  realIdx.NewVec(mol.NAtoms);
  int nat = 0;
  for (int i = 0; i != mol.NAtoms; i++) {
    realIdx(i) = nat;
    nat++;
  }

  for (int i = 0; i < nfuncs; i++) {
    std::string func = funcs[i];
    d4par(func, par, true);

    energy = 0.0;
    info =
      get_dispersion(mol, realIdx, charge, d4, par, cutoff, energy, nullptr);
    if (info != EXIT_SUCCESS) return info;

    if (check(energy, ref[i]) == EXIT_FAILURE) {
      print_fail(funcs[i], energy, ref[i]);
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

int test_param() {
  int info;

  TCutoff cutoff = TCutoff(disp2_default, 15.0);
  info = test_rational_damping(ref, cutoff);
  if (info != EXIT_SUCCESS) return info;

  // do not use cutoffs
  double huge = std::numeric_limits<double>::max();
  TCutoff no_cutoff = TCutoff(huge, huge, huge, huge);
  info = test_rational_damping(ref_no_cutoff, no_cutoff);
  if (info != EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
