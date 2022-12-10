#include <iostream>

#include <dftd_damping.h>
#include <dftd_dispersion.h>

#include "molecules.h"
#include "test_disp.h"
#include "util.h"

using namespace dftd;


int test_energy( 
  const int n,
  const char atoms[][3],
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

  // dispersion main function
  info = DFTVDW_D4(mol, par, charge, energy, nullptr);
  if (!info == EXIT_SUCCESS) return info;

  return check(energy, ref);
};

int test_disp() {
  int info;
  
  info = test_energy(16, mb16_43_01_atoms, mb16_43_01_coord, mb16_43_01_charge, mb16_43_01_ref_energy);
  if (!info == EXIT_SUCCESS) return info;

  info = test_energy(22, rost61_m1_atoms, rost61_m1_coord, rost61_m1_charge, rost61_m1_ref_energy);
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
