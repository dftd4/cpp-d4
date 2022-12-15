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
#include <dftd_geometry.h>
#include <dftd_matrix.h>
#include <dftd_ncoord.h>
#include <dftd_readxyz.h>

#include "molecules.h"
#include "test_ncoord.h"
#include "util.h"

using namespace dftd;


int test_cn(
  const int n,
  const char atoms[][4],
  const double coord[],
  const double ref_cn[]
) {
  int info;

  // assemble molecule
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (!info == EXIT_SUCCESS) return info;

  // get reference
  TVector<double> ref;
  ref.New(n);
  for (int i = 0; i != n; i++) {
    ref(i) = ref_cn[i];
  }

  // distances
  TMatrix<double> dist;      
  dist.New(n, n);
  calc_distances(mol, dist);

  // erf-CN without cutoff
  TVector<double> cn;
  TMatrix<double> dcndr; // empty because no gradient needed
  cn.New(n);
  info = get_ncoord_d4(mol, dist, 9999.9, cn, dcndr, false);
  if (!info == EXIT_SUCCESS) return info;

  // compare to ref
  for (int i = 0; i != n; i++) {
    if (check(cn(i), ref(i)) == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

int test_ncoord() {
  int info;
  
  info = test_cn(
    mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord, mb16_43_01_ref_cn
  );
  if (!info == EXIT_SUCCESS) return info;

  info = test_cn(
    rost61_m1_n, rost61_m1_atoms, rost61_m1_coord, rost61_m1_ref_cn
  );
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
