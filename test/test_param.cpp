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
#include <string>
#include <cmath>

#include <dftd_damping.h>
#include <dftd_geometry.h>
#include <dftd_parameters.h>

#include "molecules.h"
#include "test_param.h"
#include "util.h"

using namespace dftd;


int test_dftd4_energy(
  const TMolecule &mol,
  const int &charge,
  const dparam &par,
  double &energy
) {
  int info;

  // dispersion main function
  return DFTVDW_D4(mol, par, charge, energy, nullptr);
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}

int test_param() {
  int info;
  int charge;
  double energy;
  dparam par;

  // assemble molecule
  TMolecule mol;
  info = get_molecule(upu23_0a_n, upu23_0a_atoms, upu23_0a_coord, mol);
  if (!info == EXIT_SUCCESS) return info;
  charge = upu23_0a_charge;

  for (int i = 0; i < nfuncs; i++) {
    std::string func = funcs[i];
    d4par(func, par, true);

    energy = 0.0;
    info = DFTVDW_D4(mol, par, charge, energy, nullptr);
    if (!info == EXIT_SUCCESS) return info;

    if (check(energy, ref_energy[i]) == EXIT_FAILURE) {
      print_fail(funcs[i], energy, ref_energy[i]);
      printf("a1 = %.7f\n", par.a1);
      return EXIT_FAILURE;
    }

  }
 

  return EXIT_SUCCESS;
};
