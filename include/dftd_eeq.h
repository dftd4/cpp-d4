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

/**
 * Electronegativity equilibration (EEQ) model for DFT-D4
 */

#pragma once

#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_ncoord.h"

namespace multicharge {
using dftd4::NCoordErf;
using dftd4::NCoordErfEN;
using dftd4::TIVector;
using dftd4::TMatrix;
using dftd4::TMolecule;
using dftd4::TVector;

class ChargeModel {
  public:
    NCoordErf ncoord_erf; // coordination number
    // Constructor
    ChargeModel();
    // Virtual destructor
    virtual ~ChargeModel() {}

    /**
     * Get the EEQ charges for a given molecule.
     *
     * @param mol     The molecule object
     * @param dist    The distance matrix
     * @param charge  The total charge of the molecule
     * @param cutoff  The cutoff for the EEQ coordination number
     * @param q       The EEQ charges
     * @param dqdr    The derivative of the EEQ charges
     * @param lgrad   Flag for the gradient
     *
     * @return 0 if successful, 1 otherwise
     */
    int get_charges(
      const TMolecule &mol,
      const TMatrix<double> &dist,
      int charge,
      double cutoff,
      TVector<double> &q,
      TMatrix<double> &dqdr,
      bool lgrad
    );

    /**
     * Get the EEQ charges for a given molecule for the atoms specified by the
     * indices in `realIdx`.
     *
     * @param mol     The molecule object
     * @param realIdx The real atom indices (for excluding dummy atoms)
     * @param dist    The distance matrix
     * @param charge  The total charge of the molecule
     * @param cutoff  The cutoff for the EEQ coordination number
     * @param q       The EEQ charges
     * @param dqdr    The derivative of the EEQ charges
     * @param lgrad   Flag for the gradient
     *
     * @return 0 if successful, 1 otherwise
     */
    int get_charges(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      int charge,
      double cutoff,
      TVector<double> &q,
      TMatrix<double> &dqdr,
      bool lgrad
    );

    /**
     * Solves the Electronegativity Equalization (EEQ) equations for a molecule,
     * yielding atomic partial charges and, optionally, their derivatives.
     *
     * @param mol        Molecule object containing atomic information.
     * @param realIdx    Mapping of real atom indices (excludes dummy atoms).
     * @param dist       Interatomic distance matrix.
     * @param cn         Coordination numbers for each atom.
     * @param dcndr      Derivative of coordination numbers with respect to
     *                   atomic positions.
     * @param charge     Target total charge of the molecule.
     * @param qvec       Output vector of atomic charges (size: number of
     * atoms).
     * @param dqdr       Output matrix of charge derivatives with respect to
     *                   atomic positions. Only filled if `lgrad` is true.
     * @param lgrad      If true, compute charge gradients (dq/dr) in addition
     * to charges.
     * @param lverbose   If true, print diagnostic output (charges, EN values,
     * diagonals).
     *
     * @returns EXIT_SUCCESS (0) on success, or a nonzero error code if matrix
     *          inversion or setup fails.
     */
    int eeq_chrgeq(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      const int &charge,
      TVector<double> &qvec,
      TMatrix<double> &dqdr,
      bool lgrad = false,
      bool lverbose = false
    );

    /**
     * Computes the right-hand side vector (Xvec) for the EEQ charge equations.
     * Optionally computes its derivatives with respect to atomic positions.
     *
     * @param mol        Molecule object containing atomic information.
     * @param realIdx    Mapping of real atom indices (excludes dummy atoms).
     * @param charge     Target total charge of the molecule.
     * @param dist       Interatomic distance matrix.
     * @param cn         Coordination numbers for each atom.
     * @param dcndr      Derivative of the coordination number w.r.t. atom
     * positions
     * @param xvec       Output vector representing the effective driving term
     * for the EEQ system (size: number of atoms).
     * @param dxvecdr    Output vector of derivatives of xvec w.r.t. atom
     * positions. Only filled if `lgrad` is true.
     * @param qloc       Local charge (EEQ-BC specific).
     * @param dqlocdr    Derivative of the local charge (EEQ-BC specific).
     * @param cmat       Capacitance matrix.
     * @param dcmatdr    Derivative of the capacitance matrix.
     * @param lgrad      If true, compute derivatives of Xvec.
     *
     * @returns EXIT_SUCCESS (0) on success, or a nonzero error code if the
     * computation fails.
     */
    virtual int get_vrhs(
      const TMolecule &mol,
      const TIVector &realIdx,
      const int &charge,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TVector<double> &xvec,
      TMatrix<double> &dxvecdr,
      TVector<double> &qloc,
      TMatrix<double> &dqlocdr,
      TMatrix<double> &cmat,
      TMatrix<double> &dcmatdr,
      bool lgrad
    ) = 0;

    /**
     * Computes the Coulomb interaction matrix for the charge model.
     * The matrix represents the electrostatic interactions between atoms,
     * optionally including model-specific self-terms and normalization factors.
     * It forms the left-hand side (A) of the EEQ linear system: A Q = X.
     *
     * @param mol     Molecule object containing atomic information.
     * @param realIdx Mapping of real atom indices (excludes dummy atoms).
     * @param dist    Interatomic distance matrix.
     * @param cn      Coordination numbers for each atom.
     * @param qloc    Local charge (EEQ-BC specific)
     * @param cmat    Capacitance matrix (EEQ-BC specific)
     * @param Amat    Output Coulomb matrix (size: number of atoms + 1 for
     * charge constraint).
     *
     * @returns EXIT_SUCCESS (0) on success, or a nonzero error code if the
     * computation fails.
     */
    virtual int get_amat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TVector<double> &qloc,
      const TMatrix<double> &cmat,
      TMatrix<double> &Amat
    ) const = 0;

    /**
     * Computes the derivatives of the Coulomb interaction matrix with respect
     * to atomic positions. These derivatives are used to propagate the effect
     * of the charge distribution on forces or gradients in the EEQ linear
     * system.
     *
     * @param mol     Molecule object containing atomic information.
     * @param realIdx Mapping of real atom indices (excludes virtual atoms).
     * @param dist    Interatomic distance matrix.
     * @param q       Atomic charges computed from the EEQ model.
     * @param Amat    Output coulomb matrix.
     * @param dAmatdr   Output derivative matrix with respect to Cartesian
     * coordinates.
     * @param atrace  Output trace contributions for the derivative.
     * @param cn         Coordination numbers for each atom.
     * @param dcndr      Derivative of the coordination number w.r.t. atom
     * positions
     * @param qloc       Local charge (EEQ-BC specific).
     * @param dqlocdr    Derivative of the local charge (EEQ-BC specific).
     * @param cmat       Capacitance matrix.
     * @param dcmatdr    Derivative of the capacitance matrix.
     *
     * @returns EXIT_SUCCESS (0) on success, or a nonzero error code if the
     * computation fails.
     */
    virtual int get_damat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &q,
      const TMatrix<double> &Amat,
      TMatrix<double> &dAmatdr,
      TMatrix<double> &atrace,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      const TVector<double> &qloc,
      const TMatrix<double> &dqlocdr,
      const TMatrix<double> &cmat,
      const TMatrix<double> &dcmatdr
    ) const = 0;

    /**
     * This function calculates the coordination number (CN) for each atom in
     * the molecule. It also optionally computes the derivatives of CN with
     * respect to Cartesian coordinates if `lgrad` is true. The actual
     * computation is forwarded to the `get_ncoord` method of the `NCoordErf`
     * object specific to the derived charge model.
     *
     * @param mol       Molecule object containing atomic information.
     * @param realIdx   Mapping of real atom indices (excluding dummy atoms).
     * @param dist      Interatomic distance matrix.
     * @param cn        Output vector of coordination numbers for each atom.
     * @param dcndr     Output matrix of derivatives of CN with respect to
     * Cartesian coordinates.
     * @param lgrad     If true, compute derivatives of CN; otherwise, skip
     * derivative calculation.
     *
     * @returns EXIT_SUCCESS (0) on success, or a nonzero error code if the
     * computation fails.
     *
     * @throws std::runtime_error if the coordination number computation fails
     * in the derived class.
     */
    virtual int get_cn(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TVector<double> &cn,
      TMatrix<double> &dcndr,
      bool lgrad
    ) = 0;
};

// Derived class for EEQ charge model
class EEQModel : public ChargeModel {
  public:
    const double *xi;    // Element-specific electronegativity
    const double *gam;   // Element-specific chemical hardnesses
    const double *kappa; // Element-specific CN scaling constant
    const double *alp;   // Element-specific atomic radii

    // Constructs an EEQModel with element-specific parameters and
    // initializes its coordination number function.
    EEQModel();

    // Computes the right-hand side vector (xvec) for the EEQ charge equations.
    // Optionally computes its derivatives with respect to atomic positions.
    int get_vrhs(
      const TMolecule &mol,
      const TIVector &realIdx,
      const int &charge,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TVector<double> &xvec,
      TMatrix<double> &dxvecdr,
      TVector<double> &qloc,
      TMatrix<double> &dqlocdr,
      TMatrix<double> &cmat,
      TMatrix<double> &dcmatdr,
      bool lgrad
    ) override;

    // Calculate the Coulomb matrix
    int get_amat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TVector<double> &qloc,
      const TMatrix<double> &cmat,
      TMatrix<double> &Amat
    ) const override;

    // Calculate the Coulomb matrix derivatives
    int get_damat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &q,
      const TMatrix<double> &Amat,
      TMatrix<double> &dAmatdr,
      TMatrix<double> &atrace,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      const TVector<double> &qloc,
      const TMatrix<double> &dqlocdr,
      const TMatrix<double> &cmat,
      const TMatrix<double> &dcmatdr
    ) const override;

    // Calculate the coordination number, forwarding to get_ncoord
    int get_cn(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TVector<double> &cn,
      TMatrix<double> &dcndr,
      bool lgrad
    ) override;
};

// Derived class for EEQ-BC charge model
class EEQBCModel : public ChargeModel {
  private:
    static constexpr double ncoorderf_kcn = 2.0;
    static constexpr double ncoorderf_norm_exp = 0.75;
    static constexpr double ncoorderf_cutoff = 25.0;
    static constexpr double ncoorderf_f_directed = 1.0;
    static constexpr double ncoorderfen_f_directed = -1.0;
    static constexpr double ncoorderf_cn_max = -1.0;

  public:
    const double kcnrad;
    const double kbc;    // scaling factor in erf() of bond capacitance matrix
    const double cutoff; // coordination number distance cutoff
    const double cn_exp;
    const double norm_exp;
    const double *chi; // Element-specific electronegativity
    const double *eta; // Element-specific chemical hardnesses
    const double *rad; // Element-specific charge widths
    const double
      *kcnchi; // Element-specific CN scaling of the electronegativity
    const double
      *kqchi; // Element-specific local q scaling of the electronegativity
    const double
      *kqeta; // Element-specific local q scaling of the chemical hardness
    const double *cap;       // Element-specific bond capacitance
    const double *cov_radii; // Element-specific covalent radii for the CN
    const double
      *avg_cn; // Element-specific averaged coordination number over the fitset
    const double
      *rvdw; // Element-pair-specific van der Waals radii based on distance
    TVector<double> qloc; // local charges
    TMatrix<double> cmat; // capacitance matrix

    // Constructs an EEQ-BC model with element-specific parameters
    // and initializes its coordination number function.
    EEQBCModel();

    // Computes the right-hand side vector (Xvec) for the EEQ charge equations.
    // Optionally computes its derivatives with respect to atomic positions.
    int get_vrhs(
      const TMolecule &mol,
      const TIVector &realIdx,
      const int &charge,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TVector<double> &xvec,
      TMatrix<double> &dxvecdr,
      TVector<double> &qloc,
      TMatrix<double> &dqlocdr,
      TMatrix<double> &cmat,
      TMatrix<double> &dcmatdr,
      bool lgrad
    ) override;

    // Calculate the Coulomb matrix
    int get_amat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TVector<double> &qloc,
      const TMatrix<double> &cmat,
      TMatrix<double> &Amat
    ) const override;

    // Calculate the Coulomb matrix derivatives
    int get_damat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &q,
      const TMatrix<double> &Amat,
      TMatrix<double> &dAmatdr,
      TMatrix<double> &atrace,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      const TVector<double> &qloc,
      const TMatrix<double> &dqlocdr,
      const TMatrix<double> &cmat,
      const TMatrix<double> &dcmatdr
    ) const override;

    // Calculate the coordination number, forwarding to get_ncoord
    int get_cn(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TVector<double> &cn,
      TMatrix<double> &dcndr,
      bool lgrad
    ) override;

    // Get purely geometry-dependent local charges
    int get_qloc(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const double q_tot,
      TVector<double> &qloc,
      TMatrix<double> &dqlocdr,
      bool lgrad
    );

    // Get the capacitance for bond between atoms i and j
    int get_cpair(int iat, int jat, double &dist_ij, double &c_ij) const;

    // Get the capacitance for bond between atoms i and j
    int get_dcpair(
      double dist_ij,   // distance between atom j to atom i
      double rvdw_ijat, // pairwise van der Waals radii for ij atom types
      double cap_ij, // product of bond capacitances for atom types of i and j
      TVector<double> &vec,    // Vector from j to i
      TVector<double> &dcdr_ij // Out: Capacitance for bond ij
    ) const;

    // Get the capacitance matrix
    int get_cmat(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TMatrix<double> &cmat
    );

    // Get the capacitance matrix
    int get_dcmatdr(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TMatrix<double> &cmat
    );

    // Get the right-hand side (electronegativity) of the set of linear
    // equations
    int get_xvec(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      TMatrix<double> &cmat,
      int charge,
      TVector<double> &qloc,
      TVector<double> &xvec
    );

    // Get the derivative of the right-hand side of the set of linear equations
    int get_xvec_derivs(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TMatrix<double> &cmat,
      int charge,
      TVector<double> &xvec,
      TMatrix<double> &dxvecdr,
      TVector<double> &qloc,
      TMatrix<double> &dqlocdr,
      TMatrix<double> &dcmatdr
    );
};
} // namespace multicharge
