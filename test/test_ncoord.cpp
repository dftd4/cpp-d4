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

using namespace dftd4;

int test_cn(
  int n,
  const char atoms[][3],
  const double coord[],
  const double ref_cn[]
) {
  int info;

  // assemble molecule
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (info != EXIT_SUCCESS) return info;

  // get reference
  TVector<double> ref;
  ref.New(n);
  for (int i = 0; i != n; i++) {
    ref(i) = ref_cn[i];
  }

  // dummy for masking of possible ghost atoms
  dftd4::TIVector realIdx;
  realIdx.NewVec(mol.NAtoms);
  int nat = 0;
  for (int i = 0; i != mol.NAtoms; i++) {
    realIdx(i) = nat++;
  }

  // distances
  TMatrix<double> dist;
  dist.New(n, n);
  calc_distances(mol, realIdx, dist);

  // erf-CN without cutoff
  TVector<double> cn;
  TMatrix<double> dcndr;
  NCoordErfD4 ncoord_erf(7.5, 1.0, 9999.9);
  info = ncoord_erf.get_ncoord(mol, realIdx, dist, cn, dcndr, false);
  if (info != EXIT_SUCCESS) return info;

  // compare to ref
  for (int i = 0; i != n; i++) {
    if (check(cn(i), ref(i)) == EXIT_FAILURE) {
      print_fail("CN_D4", cn(i), ref(i));
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

int test_numgrad_d4(int n, const char atoms[][3], const double coord[]) {
  // assemble molecule
  int info;
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (info != EXIT_SUCCESS) return info;
  double step{1.0e-6};
  double thr{1e-8};

  TMatrix<double> dist;
  dist.NewMat(mol.NAtoms, mol.NAtoms);
  TMatrix<double> num_dcndr; // numerical gradient of the coordination number
  num_dcndr.NewMat(mol.NAtoms, 3 * mol.NAtoms);
  TMatrix<double>
    analytic_dcndr; // analytical gradient of the coordination number
  analytic_dcndr.NewMat(mol.NAtoms, 3 * mol.NAtoms);
  NCoordErfD4 ncoord_erf_d4;
  TVector<double> cn, cn_r, cn_l;
  TMatrix<double> dcndr_dum;

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
  ncoord_erf_d4.get_ncoord(mol, realIdx, dist, cn, analytic_dcndr, true);

  // check if analytical gradient is antisymmetric
  for (int c = 0; c < 3; c++) {
    for (int i = 0; i < mol.NAtoms; i++) {
      for (int k = 0; k < i; k++) {
        if (abs(analytic_dcndr(k, 3 * i + c) + analytic_dcndr(i, 3 * k + c)) >
            1.0e-9) {
          print_fail(
            "Analytical CN-gradient is not antisymmetric for NCoordErfD4",
            analytic_dcndr(k, 3 * i + c) + analytic_dcndr(i, 3 * k + c),
            0.0
          );
          return EXIT_FAILURE;
        }
      }
    }
  }

  // numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      mol.CC(i, c) += step;
      calc_distances(mol, realIdx, dist);
      ncoord_erf_d4.get_ncoord(mol, realIdx, dist, cn_r, dcndr_dum, false);

      mol.CC(i, c) = mol.CC(i, c) - 2 * step;
      calc_distances(mol, realIdx, dist);
      ncoord_erf_d4.get_ncoord(mol, realIdx, dist, cn_l, dcndr_dum, false);

      mol.CC(i, c) = mol.CC(i, c) + step;
      for (int j = 0; j < mol.NAtoms; j++) {
        // numerical CN gradient: dCN(j)/ dr(i)^c with c = x,y, or z
        num_dcndr(j, 3 * i + c) = 0.5 * (cn_r(j) - cn_l(j)) / step;
      }
    }
  }

  // compare against numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      for (int j = 0; j < mol.NAtoms; j++) {
        if (check(analytic_dcndr(j, 3 * i + c), num_dcndr(j, 3 * i + c), thr) !=
            EXIT_SUCCESS) {
          print_fail(
            "Gradient mismatch for NCoordErfD4 dcndr",
            analytic_dcndr(j, 3 * i + c),
            num_dcndr(j, 3 * i + c)
          );
          return EXIT_FAILURE;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

int test_numgrad(int n, const char atoms[][3], const double coord[]) {
  // assemble molecule
  int info;
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (info != EXIT_SUCCESS) return info;
  double step{1.0e-6};
  double thr{1e-8};

  TMatrix<double> dist;
  dist.NewMat(mol.NAtoms, mol.NAtoms);
  TMatrix<double> num_dcndr; // numerical gradient of the coordination number
  num_dcndr.NewMat(mol.NAtoms, 3 * mol.NAtoms);
  TMatrix<double>
    analytic_dcndr; // analytical gradient of the coordination number
  analytic_dcndr.NewMat(mol.NAtoms, 3 * mol.NAtoms);
  NCoordErf ncoord_erf;
  TVector<double> cn, cn_r, cn_l;
  TMatrix<double> dcndr_dum;

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
  ncoord_erf.get_ncoord(mol, realIdx, dist, cn, analytic_dcndr, true);

  // numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      mol.CC(i, c) += step;
      calc_distances(mol, realIdx, dist);
      ncoord_erf.get_ncoord(mol, realIdx, dist, cn_r, dcndr_dum, false);

      mol.CC(i, c) = mol.CC(i, c) - 2 * step;
      calc_distances(mol, realIdx, dist);
      ncoord_erf.get_ncoord(mol, realIdx, dist, cn_l, dcndr_dum, false);

      mol.CC(i, c) = mol.CC(i, c) + step;
      for (int j = 0; j < mol.NAtoms; j++) {
        // numerical CN gradient: dCN(j)/ dr(i)^c with c = x,y, or z
        num_dcndr(j, 3 * i + c) = 0.5 * (cn_r(j) - cn_l(j)) / step;
      }
    }
  }

  // compare against numerical gradient
  for (int i = 0; i < mol.NAtoms; i++) {
    for (int c = 0; c < 3; c++) {
      for (int j = 0; j < mol.NAtoms; j++) {
        if (check(analytic_dcndr(j, 3 * i + c), num_dcndr(j, 3 * i + c), thr) !=
            EXIT_SUCCESS) {
          print_fail(
            "Gradient mismatch for NCoordErf dcndr",
            analytic_dcndr(j, 3 * i + c),
            num_dcndr(j, 3 * i + c)
          );
          return EXIT_FAILURE;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

int test_ncoord() {
  int info;

  info = test_cn(water_n, water_atoms, water_coord, water_ref_cn);
  if (info != EXIT_SUCCESS) return info;

  info = test_cn(
    mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord, mb16_43_01_ref_cn
  );
  if (info != EXIT_SUCCESS) return info;

  info =
    test_cn(rost61_m1_n, rost61_m1_atoms, rost61_m1_coord, rost61_m1_ref_cn);
  if (info != EXIT_SUCCESS) return info;

  info = test_numgrad_d4(mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord);
  if (info != EXIT_SUCCESS) return info;

  info = test_numgrad(mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord);
  if (info != EXIT_SUCCESS) return info;

  info = test_numgrad_d4(water_n, water_atoms, water_coord);
  if (info != EXIT_SUCCESS) return info;

  info = test_numgrad(water_n, water_atoms, water_coord);
  if (info != EXIT_SUCCESS) return info;
  return EXIT_SUCCESS;
}
