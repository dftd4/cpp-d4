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
  using dftd4::TIVector;
  using dftd4::TVector;
  using dftd4::TMatrix;
  using dftd4::TMolecule;
  using dftd4::NCoordErf;
  using dftd4::NCoordErfEN;

class ChargeModel
{
public:
  NCoordErf ncoord_erf;     // coordination number
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
  
  virtual int get_vrhs(
    const TMolecule &mol,
    const TIVector &realIdx,
    const int &charge,
    const TMatrix<double> &dist,
    const TVector<double> &cn,
    const TMatrix<double> &dcndr,
    TVector<double> &Xvec,
    TVector<double> &dXvec,
    bool lgrad
  ) = 0;

  // Calculate the Coulomb matrix
  virtual int get_amat_0d(
    const TMolecule &mol,
    const TIVector &realIdx,
    const TMatrix<double> &dist,
    const TVector<double> &cn,
    const TMatrix<double> &dcndr,
    TMatrix<double> &Amat
  ) const = 0;
  
  // Calculate the Coulomb matrix derivatives
  virtual int get_damat_0d(
    const TMolecule &mol,
    const TIVector &realIdx,
    const TMatrix<double> &dist,
    const TVector<double> &q,
    const TMatrix<double> &Amat,
    TMatrix<double> &dAmat,
    TMatrix<double> &atrace
  ) const = 0;

  // Calculate the coordination number, forwarding to get_ncoord
  virtual NCoordErf& get_cn(
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
    const double* xi;     // Element-specific electronegativity 
    const double* gam;    // Element-specific chemical hardnesses
    const double* kappa;  // Element-specific CN scaling constant
    const double* alp;    // Element-specific atomic radii
    
    //
    EEQModel();
  
    int get_vrhs(
      const TMolecule &mol,
      const TIVector &realIdx,
      const int &charge,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TVector<double> &Xvec,
      TVector<double> &dXvec,
      bool lgrad
    ) override;
  
    // Calculate the Coulomb matrix
    int get_amat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TMatrix<double> &Amat
    ) const override;
    
    // Calculate the Coulomb matrix derivatives
    int get_damat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &q,
      const TMatrix<double> &Amat,
      TMatrix<double> &dAmat,
      TMatrix<double> &atrace
    ) const override;

    // Calculate the coordination number, forwarding to get_ncoord
    NCoordErf& get_cn(
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
    const double kbc;         // scaling factor in erf() of bond capacitance matrix
    const double cutoff;      // coordination number distance cutoff
    const double cn_exp;
    const double norm_exp;
    const double* chi;        // Element-specific electronegativity
    const double* eta;        // Element-specific chemical hardnesses 
    const double* rad;        // Element-specific charge widths
    const double* kcnchi;     // Element-specific CN scaling of the electronegativity
    const double* kqchi;      // Element-specific local q scaling of the electronegativity
    const double* kqeta;      // Element-specific local q scaling of the chemical hardness
    const double* cap;        // Element-specific bond capacitance
    const double* cov_radii;  // Element-specific covalent radii for the CN 
    const double* avg_cn;     // Element-specific averaged coordination number over the fitset
    const double* rvdw;       // Element-pair-specific van der Waals radii based on distance
    TVector<double> qloc;     // local charges
    TMatrix<double> cmat;     // capacitance matrix

    //
    EEQBCModel();
  
    int get_vrhs(
      const TMolecule &mol,
      const TIVector &realIdx,
      const int &charge,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TVector<double> &Xvec,
      TVector<double> &dXvec,
      bool lgrad
    ) override;
  
    // Calculate the Coulomb matrix
    int get_amat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TMatrix<double> &Amat
    ) const override;
    
    // Calculate the Coulomb matrix derivatives
    int get_damat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      const TVector<double> &q,
      const TMatrix<double> &Amat,
      TMatrix<double> &dAmat,
      TMatrix<double> &atrace
    ) const override;

    // Calculate the coordination number, forwarding to get_ncoord
    NCoordErf& get_cn(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TVector<double> &cn,
      TMatrix<double> &dcndr,
      bool lgrad
    ) override;

    // Get purely geometry-dependent local charges
    int get_qloc(
      const TMolecule&,
      const TIVector &,
      const TMatrix<double>&,
      const double,
      TVector<double>&
    );

    // Get the capacitance for bond between atoms i and j
    int get_cpair(
      int iat,
      int jat,
      double &dist_ij,
      double &c_ij
    ) const;

    // Get the capacitance matrix
    int get_cmat(
      const TMolecule&,
      const TIVector &realIdx,
      const TMatrix<double>&,
      TMatrix<double>&
    );

    // Get the right-hand side (electronegativity) of the set of linear equations
    int get_xvec(
      const TMolecule&,
      const TIVector&,
      const TMatrix<double>&,
      const TVector<double> &cn,
      const TMatrix<double> &dcndr,
      TMatrix<double>&,
      int,
      TVector<double>&
    );

    // numerical gradient of partial charges w.r.t. atom positions
    int num_grad_dqdr(
      TMolecule&,  // molecular geometry
      const TIVector&,
      int,
      TMatrix<double>& // numerical gradient
    );

};
} // namespace multicharge
