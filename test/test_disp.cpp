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

#include "molecules.h"
#include "test_disp.h"
#include "util.h"

using namespace dftd;


int test_energy( 
  const int n,
  const char atoms[][4],
  const double coord[],
  const int charge,
  const double ref
) {
  int info = 0;
  bool lmbd = true;
  double energy = 0.0;

  // BP86 parameters
  dparam par;
  d4par("bp86", par, lmbd);

  // assemble molecule
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (!info == EXIT_SUCCESS) return info;

  TCutoff cutoff;

  // dispersion main function
  info = get_dispersion(mol, par, charge, cutoff, energy, nullptr);
  if (!info == EXIT_SUCCESS) return info;

  if (check(energy, ref) == EXIT_FAILURE) {
    print_fail("BP86-D4-ATM", energy, ref);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int test_disp() {
  int info;

  info = test_energy(water_n, water_atoms, water_coord, water_charge, water_ref_energy);
  if (!info == EXIT_SUCCESS) return info;
  
  info = test_energy(16, mb16_43_01_atoms, mb16_43_01_coord, mb16_43_01_charge, mb16_43_01_ref_energy);
  if (!info == EXIT_SUCCESS) return info;

  info = test_energy(22, rost61_m1_atoms, rost61_m1_coord, rost61_m1_charge, rost61_m1_ref_energy);
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
