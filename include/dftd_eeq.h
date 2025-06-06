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

namespace multicharge {
  using dftd4::TIVector;
  using dftd4::TVector;
  using dftd4::TMatrix;
  using dftd4::TMolecule;

class ChargeModel
{
public:
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
    const int &charge,
    const TVector<double> &cn,
    TVector<double> &qvec,
    TMatrix<double> &dcndr,
    TMatrix<double> &dqdr,
    bool lgrad = false,
    bool lverbose = false
  );
  
  virtual int get_vrhs(
    const TMolecule &mol,
    const TIVector &realIdx,
    const int &charge,
    const TVector<double> &cn,
    TVector<double> &Xvec,
    TVector<double> &dXvec,
    bool lgrad
  ) const = 0;

  // Calculate the Coulomb matrix
  virtual int get_amat_0d(
    const TMolecule &mol,
    const TIVector &realIdx,
    const TMatrix<double> &dist,
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
  
};

// Derived class for EEQ charge model
class EEQModel : public ChargeModel {
  public:
    const double* xi;
    const double* gam;
    const double* kappa;
    const double* alp;
    
    //
    EEQModel();
  
    int get_vrhs(
      const TMolecule &mol,
      const TIVector &realIdx,
      const int &charge,
      const TVector<double> &cn,
      TVector<double> &Xvec,
      TVector<double> &dXvec,
      bool lgrad
    ) const override;
  
    // Calculate the Coulomb matrix
    int get_amat_0d(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
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
};

} // namespace multicharge
