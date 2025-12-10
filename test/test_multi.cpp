
#include "test_multi.h"
#include "dftd_eeq.h"
#include "dftd_multicharge_param.h"
#include "molecules.h"
#include "util.h"

using namespace dftd4;

int test_multi_param() {
  // Test EEQ-BC model parameters
  // Since they are passed via the constructor initializer list, the sorting has
  // to be the same
  //  as in the EEQ-BC parameter header file
  // Chi parameter (currently first in constructor initializer list)
  if (check(
        multicharge_param::eeqbc::eeqbc_chi[1], param_eeqbc_chi_ref[0], 1.0E-9
      ) == EXIT_FAILURE) {
    print_fail(
      "Multicharge: Param",
      multicharge_param::eeqbc::eeqbc_chi[1],
      param_eeqbc_chi_ref[0]
    );
    return EXIT_FAILURE;
  }
  if (check(
        multicharge_param::eeqbc::eeqbc_chi[103], param_eeqbc_chi_ref[1], 1.0E-9
      ) == EXIT_FAILURE) {
    print_fail(
      "Multicharge: Param",
      multicharge_param::eeqbc::eeqbc_chi[103],
      param_eeqbc_chi_ref[1]
    );
    return EXIT_FAILURE;
  }
  // rvdw parameter (currently last in constructor initializer list)
  // Also check correct indexing
  int iat = 1;
  int jat = 1;
  int ij_at = (((iat - 1) * iat) / 2 + jat) -
              1; // calculate index for half-vectorized matrix rvdw
  if (check(
        multicharge_param::eeqbc::eeqbc_rvdw[ij_at],
        param_eeqbc_rvdw_ref[0],
        1.0E-4
      ) == EXIT_FAILURE) {
    print_fail(
      "Multicharge: Param",
      multicharge_param::eeqbc::eeqbc_rvdw[ij_at],
      param_eeqbc_rvdw_ref[0]
    );
    return EXIT_FAILURE;
  }
  iat = 103;
  jat = 102;
  ij_at = (((iat - 1) * iat) / 2 + jat) -
          1; // calculate index for half-vectorized matrix rvdw
  if (check(
        multicharge_param::eeqbc::eeqbc_rvdw[ij_at],
        param_eeqbc_rvdw_ref[1],
        1.0E-4
      ) == EXIT_FAILURE) {
    print_fail(
      "Multicharge: Param",
      multicharge_param::eeqbc::eeqbc_rvdw[ij_at],
      param_eeqbc_rvdw_ref[1]
    );
    return EXIT_FAILURE;
  }
  iat = 103;
  jat = 103;
  ij_at = (((iat - 1) * iat) / 2 + jat) -
          1; // calculate index for half-vectorized matrix rvdw
  if (check(
        multicharge_param::eeqbc::eeqbc_rvdw[ij_at],
        param_eeqbc_rvdw_ref[2],
        1.0E-4
      ) == EXIT_FAILURE) {
    print_fail(
      "Multicharge: Param",
      multicharge_param::eeqbc::eeqbc_rvdw[ij_at],
      param_eeqbc_rvdw_ref[2]
    );
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// Test member functions of ChargeModel and derived classes
int test_multi_functions() {
  int info;
  // assemble molecule
  TMolecule mol;
  info = get_molecule(mb16_43_01_n, mb16_43_01_atoms, mb16_43_01_coord, mol);
  if (info != EXIT_SUCCESS) {
    printf("Multicharge: Functions, Failed to set up molecule.");
    fflush(stdout);
    return info;
  }

  // Test EEQ-BC member functions
  multicharge::EEQBCModel eeqbc_model;
  TIVector realIdx;
  TMatrix<double> dist;
  dftd4::initializeRealIdx(mb16_43_01_n, realIdx);
  dist.NewMatrix(mb16_43_01_n, mb16_43_01_n);
  dftd4::calc_distances(mol, realIdx, dist);
  TVector<double> cn;
  TMatrix<double> dcndr;

  // Calculate partial charges
  TVector<double> q;
  TMatrix<double> dqdr;
  q.NewVector(mb16_43_01_n);
  dqdr.NewMatrix(3 * mb16_43_01_n, mb16_43_01_n);

  eeqbc_model.get_cn(mol, realIdx, dist, cn, dcndr, true);
  info = eeqbc_model.eeq_chrgeq(
    mol, realIdx, dist, cn, dcndr, mb16_43_01_charge, q, dqdr, true, false
  );
  if (info != EXIT_SUCCESS) {
    printf("Multicharge: Functions, Failed to calculate charges.");
    fflush(stdout);
    return info;
  }

  // Check against multicharge reference calculation
  for (int i = 0; i < mol.NAtoms; i++) {
    if (check(q(i), qvec_reference[i], 1.0E-8) == EXIT_FAILURE) {
      print_fail(
        "Multicharge: Functions, Partial charge differs from reference.",
        q(i),
        qvec_reference[i]
      );
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

// Test member functions of ChargeModel and derived classes
int test_eeqbc_amalgam() {
  int info;
  // assemble molecule
  TMolecule mol;
  info = get_molecule(amalgam_n, amalgam_atoms, amalgam_coord, mol);
  if (info != EXIT_SUCCESS) {
    printf("Multicharge: Functions, Failed to set up molecule.");
    fflush(stdout);
    return info;
  }

  // Test EEQ-BC member functions
  multicharge::EEQBCModel eeqbc_model;
  TIVector realIdx;
  TMatrix<double> dist;
  dftd4::initializeRealIdx(amalgam_n, realIdx);
  dist.NewMatrix(amalgam_n, amalgam_n);
  dftd4::calc_distances(mol, realIdx, dist);
  TVector<double> cn;
  TMatrix<double> dcndr;

  // Calculate partial charges
  TVector<double> q;
  TMatrix<double> dqdr;
  q.NewVector(amalgam_n);
  dqdr.NewMatrix(3 * amalgam_n, amalgam_n);

  eeqbc_model.get_cn(mol, realIdx, dist, cn, dcndr, true);

  info = eeqbc_model.eeq_chrgeq(
    mol, realIdx, dist, cn, dcndr, amalgam_charge, q, dqdr, false, false
  );
  if (info != EXIT_SUCCESS) {
    printf("Multicharge: Failed to calculate charges.");
    fflush(stdout);
    return info;
  }

  // Check against multicharge reference calculation
  for (int i = 0; i < mol.NAtoms; i++) {
    if (check(q(i), qvec_amalgam_reference[i], 1.0E-8) == EXIT_FAILURE) {
      print_fail(
        "Multicharge: Functions, Partial charge differs from reference for "
        "amalgam.",
        q(i),
        qvec_amalgam_reference[i]
      );
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

// Test for multicharge models
int test_multi() {
  int info;

  // Test the charge model parameters
  info = test_multi_param();
  if (info != EXIT_SUCCESS) return info;

  // Test member functions of the ChargeModel and derived classes
  info = test_multi_functions();
  if (info != EXIT_SUCCESS) return info;

  // Test amalgam structure
  info = test_eeqbc_amalgam();
  if (info != EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
}
