#include <iostream>

#include <dftd_eeq.h>
#include <dftd_damping.h>
#include <dftd_ncoord.h>

#include "molecules.h"
#include "test_disp.h"
#include "util.h"

using namespace dftd;


int test_energy( 
  const int n,
  const char atoms[][3],
  const double coord[],
  const double ref
) {
  int info = 0;
  int ndim = 0;
  int charge = 0;

  double wf = 6.0, g_a = 3.0, g_c = 2.0;
  double energy = 0.0;

  bool lmbd = true;
  bool lgrad = false;

  // BP86 parameters
  dparam par;
  d4par("bp86", par, lmbd);

  // assemble molecule
  TMolecule mol;
  info = get_molecule(n, atoms, coord, mol);
  if (!info == EXIT_SUCCESS) return info;

  // distances
  TMatrix<double> dist;      
  dist.New(n, n);
  calc_distances(mol, dist);

  // calculation dimension of D4
  ndim = d4dim(mol);

  double es = 0.0;           // electrostatic energy
  TVector<double> covcn;     // D4 coordination number
  TVector<double> cn;        // EEQ cordination number
  TVector<double> q;         // partial charges from EEQ model
  TVector<double> gweights;  // Gaussian weights for C6 interpolation
  TMatrix<double> c6ref;     // reference C6 coefficients
  TMatrix<double> numg;      // derivative of dispersion energy

  // get memory
  c6ref.New(ndim, ndim);
  covcn.New(mol.NAtoms);
  cn.New(n);
  gweights.New(ndim);
  q.New(mol.NAtoms + 1);


  TMatrix<double> dcndr;     // derivative of erf-CN
  TMatrix<double> dcovcndr;  // derivative of covalent D4
  TMatrix<double> dqdr;      // derivative of partial charges
  TMatrix<double> ges;       // derivative of electrostatic energy
  TMatrix<double> gradient;  // derivative of dispersion energy
  if (lgrad) {
    dcndr.New(mol.NAtoms, 3 * mol.NAtoms);
    dcovcndr.New(mol.NAtoms, 3 * mol.NAtoms);
    dqdr.New(mol.NAtoms + 1, 3 * mol.NAtoms);
    ges.New(mol.NAtoms, 3);
    gradient.New(mol.NAtoms, 3);
  } 
  
  // get the EEQ coordination number
  info = get_ncoord_erf(mol, dist, cn, dcndr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // calculate partial charges by EEQ model
  info = eeq_chrgeq(
    mol, charge, dist, cn, q, es, dcndr, dqdr, ges, lgrad, false, false
  );
  if (!info == EXIT_SUCCESS) return info;

  // get the D4 coordination number
  info = get_ncoord_d4(mol, dist, covcn, dcovcndr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // D4 weights and c6 references
  info = d4(mol, ndim, wf, g_a, g_c, covcn, gweights, c6ref);
  if (!info == EXIT_SUCCESS) return info;

  // D4 energy
  info = edisp(mol, dist, ndim, q, par, g_a, g_c, gweights, c6ref, lmbd, energy);
  if (!info == EXIT_SUCCESS) return info;

  return check(energy, ref);
};

int test_disp() {
  int info;
  
  info = test_energy(16, mb16_43_01_atoms, mb16_43_01_coord, mb16_43_01_ref_energy);
  if (!info == EXIT_SUCCESS) return info;

  info = test_energy(22, rost61_m1_atoms, rost61_m1_coord, rost61_m1_ref_energy);
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
