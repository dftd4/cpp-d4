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

/*
 * Definition of the coordination numbers used in DFT-D4 and the
 * electronegativity equilibration (EEQ) model
 */

#pragma once

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd4 {


class NCoordBase
{
  public:
    TVector<double> cn;
    TMatrix<double> dcndr;
    static const double rad[];
    double kcn;  // Steepness of counting function 
    double norm_exp;
    double cutoff;
    // Get the coordination number
    int get_ncoord( // with ghost atoms
      const TMolecule&,
      const TMatrix<double>&,
      const double,
      TVector<double>&,
      TMatrix<double>&,
      bool);
    int get_ncoord(  // without ghost atoms
      const TMolecule&,
      const TIVector&,
      const TMatrix<double>&,
      bool);
    // Calculate the coordination number using the virtual counting function
    int ncoord_base(
      const TMolecule&,
      const TIVector&,
      const TMatrix<double>&);
    // Calculate the derivative of the coordination number
    int dncoord_base(
      const TMolecule&,
      const TIVector&,
      const TMatrix<double>&);
    // Get the DFT-D4 coordination number
    int get_ncoord_d4(
      const TMolecule&,
      const TMatrix<double>&,
      bool);
    int get_ncoord_d4(
      const TMolecule&,
      const TIVector&,
      const TMatrix<double>&,
      bool);
    int ncoord_d4(
      const TMolecule&,
      const TIVector&,
      const TMatrix<double>&);
    int dncoord_d4(
      const TMolecule&,
      const TIVector&,
      const TMatrix<double>&);

    // Counting function
    virtual double count_fct(double) const = 0;
    // Derivative of the counting function
    virtual double dcount_fct(double) const = 0;
    // Constructor
    NCoordBase(double optional_kcn = 7.5, double optional_norm_exp = 1.0, double optional_cutoff = 25.0);
    // Virtual destructor
    virtual ~NCoordBase() {
      cn.DelVec();
      dcndr.DelMat();
    }
};

class NCoordErf : public NCoordBase {
  public:
    // erf() based counting function
    double count_fct(double) const override;
    // derivative of the erf() based counting function
    double dcount_fct(double) const override;
    // Constructor
    NCoordErf(double optional_kcn = 7.5, double optional_norm_exp = 1.0, double optional_cutoff = 25.0)
    : NCoordBase(optional_kcn, optional_norm_exp, optional_cutoff){}
    // Use default destructor; base class handles cleanup
    ~NCoordErf() override = default;
};

/**
 * Calculate all distance pairs and store in matrix.
 *
 * @param mol Molecule object.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param dist Distance matrix (inout).
 * @return Exit status.
 */
extern int calc_distances(
  const TMolecule &mol,
  const TIVector &realIdx,
  TMatrix<double> &dist
);

/**
 * Initialize real indices to all atoms in the molecule.
 *
 * @param nat Number of atoms in the molecule.
 * @param realIdx Vector to store the real indices.
 * @return void
 */
extern void initializeRealIdx(int nat, TVector<int> &realIdx);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * Wrapper for error function coordination number.
 *
 * @param mol Molecule object.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param dist Distance matrix.
 * @param cn Vector of coordination numbers.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param dcndr Derivative of coordination number.
 * @param lgrad Flag for gradient computation.
 * @return Exit status.
 */
extern int get_ncoord_erf(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad = false
);

/**
 * Wrapper for error function coordination number.
 *
 * @param mol Molecule object.
 * @param dist Distance matrix.
 * @param cn Vector of coordination numbers.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param dcndr Derivative of coordination number.
 * @param lgrad Flag for gradient computation.
 * @return Exit status.
 */
extern int get_ncoord_erf(
  const TMolecule &mol,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad = false
);

/**
 * Calculate error function coordination number.
 *
 * @param mol Molecule object.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param dist Distance matrix.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param cn Vector of coordination numbers.
 * @return Exit status.
 */
extern int ncoord_erf(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn
);

/**
 * Calculate error function coordination number and derivative
 * w.r.t. nuclear coordinates.
 *
 * @param mol Molecule object.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param dist Distance matrix.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param cn Vector of coordination numbers.
 * @param dcndr Derivative of coordination number.
 * @return Exit status.
 */
extern int dncoord_erf(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr
);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * Wrapper for error function coordination number for DFT-D4.
 *
 * @param mol Molecule object.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param dist Distance matrix.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param cn Vector of coordination numbers.
 * @return Exit status.
 */
extern int get_ncoord_d4(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad = false
);

/**
 * Wrapper for error function coordination number for DFT-D4.
 *
 * @param mol Molecule object.
 * @param dist Distance matrix.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param cn Vector of coordination numbers.
 * @return Exit status.
 */
extern int get_ncoord_d4(
  const TMolecule &mol,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad = false
);

/**
 * Calculate covalent coordination number for DFT-D4.
 *
 * @param mol Molecule object.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param dist Distance matrix.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param cn Vector of coordination numbers.
 * @return Exit status.
 */
extern int ncoord_d4(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn
);

/**
 * Calculate covalent coordination number and derivative
 * w.r.t. nuclear coordinates
 *
 * @param mol Molecule object.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param dist Distance matrix.
 * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
 * @param cn Vector of coordination numbers.
 * @param dcndr Derivative of coordination number.
 * @return Exit status.
 */
extern int dncoord_d4(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  double cutoff,
  TVector<double> &cn,
  TMatrix<double> &dcndr
);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * Error function counting function for coordination number contributions.
 *
 * @param k Steepness of the counting function.
 * @param rr Current distance over cutoff radius.
 * @return Value of the counting function.
 */
extern double erf_count(double k, double rr);

/**
 * Derivative of the counting function w.r.t. the distance.
 *
 * @param k Steepness of the counting function.
 * @param rr TCutoff radius.
 * @return Derivative of the counting function.
 */
extern double derf_count(double k, double rr);

/**
 * TCutoff function for large coordination numbers
 *
 * @param cn_max Maximum CN (not strictly obeyed).
 * @param cn On input coordination number, on output modified CN.
 * @param dcndr On input derivative of CN w.r.t. cartesian coordinates,
 * on output derivative of modified CN.
 * @param lgrad Flag for gradient.
 */
extern int cut_coordination_number(
  double cn_max,
  TVector<double> &cn,
  TMatrix<double> &dcndr,
  bool lgrad
);

extern inline double log_cn_cut(double cn_max, double cn);

extern inline double dlog_cn_cut(double cn_max, double cn);

} // namespace dftd4
