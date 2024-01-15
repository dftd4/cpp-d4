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
#include <dftd_model.h>

#include "molecules.h"
#include "test_disp2.h"
#include "util.h"

using namespace dftd4;

int test_energy2(
  const int n,
  const char atoms[][4],
  const double coord[],
  const int charge,
  const double ref
) {
  int info{0};
  bool latm{true};
  double energy{0.0};

  // BP86 parameters
  dparam par;
  d4par("bp86", par, latm);
  par.s9 = 0.0;

  // assemble molecule
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (!info == EXIT_SUCCESS) return info;

  TCutoff cutoff;
  TD4Model d4;

  // dispersion main function
  info = get_dispersion(mol, charge, d4, par, cutoff, energy, nullptr);
  if (!info == EXIT_SUCCESS) return info;

  if (check(energy, ref, 1e-8) == EXIT_FAILURE) {
    print_fail("BP86-D4", energy, ref);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int test_disp2() {
  int info;

  info = test_energy2(
    water_n, water_atoms, water_coord, water_charge, -2.3162150393943E-04
  );
  if (!info == EXIT_SUCCESS) return info;

  info = test_energy2(
    mb16_43_01_n,
    mb16_43_01_atoms,
    mb16_43_01_coord,
    mb16_43_01_charge,
    -2.5912431304617E-02
  );
  if (!info == EXIT_SUCCESS) return info;

  info = test_energy2(
    actinides_n,
    actinides_atoms,
    actinides_coord,
    actinides_charge,
    -2.760462533236626E-01
  );
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
