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
#include "dftd_multicharge_param.h"

namespace dftd4 {

class NCoordBase
{
  public:
    TVector<double> cn;  // coordination number
    TMatrix<double> dcndr;  // derivative of the coordination number
    static const double rad[];  // cov. radii default: D3 covalent radii from dftd_ncoord.cpp
    double kcn;  // Steepness of counting function 
    double norm_exp;  // exponent of the normalizing-factor in the counting function
    double cutoff;  // Coordination number cutoff distance
    double f_directed;  // directed factor for EN scaled coordination number
    double cn_max;  // (Soft-)maximum value for the coordination number
    const double* rcov;  // covalent radii for CN
    // Get the coordination number
    /**
    * Wrapper for coordination number calculation.
    *
    * @param mol Molecule object.
    * @param dist Distance matrix.
    * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
    * @param lgrad Flag for gradient computation.
    * @return Exit status.
    */
    int get_ncoord( // without ghost atom indices
      const TMolecule &mol,
      const TMatrix<double> &dist,
      const double cutoff,
      bool lgrad);
          /**
    * Wrapper for coordination number calculation.
    *
    * @param mol Molecule object.
    * @param realIdx List for real atoms excluding ghost/non atoms.
    * @param dist Distance matrix.
    * @param cutoff Real-space cutoff (default: @see {dftd_cutoff.h}).
    * @param cn Vector of coordination numbers.
    * @return Exit status.
    */
    int get_ncoord(  // with ghost atoms
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      bool lgrad);
    /**
    * Calculate the coordination number.
    *
    * @param mol Molecule object.
    * @param realIdx List for real atoms excluding ghost/non atoms.
    * @param dist Distance matrix.
    * @return Exit status.
    */
    int ncoord_base(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist);
    /**
    * Calculate error function coordination number and derivative
    * w.r.t. nuclear coordinates.
    *
    * @param mol Molecule object.
    * @param realIdx List for real atoms excluding ghost/non atoms.
    * @param dist Distance matrix.
    * @return Exit status.
    */
    int dr_ncoord_base(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist);
    /**
    * Electronegativity factor.
    *
    * @param i atom index
    * @param j atom index
    * @return Value of the electronegativity factor.
    */
    virtual double get_en_factor(int i, int j) const;
    // Calculate the element pair-specific covalent radii
    virtual double get_rcov(int, int) const;
    /**
    * base class function for coordination number contributions.
    *
    * @param r distance between the two atoms 
    * @param rc summed up covalent radii of the two atoms
    * @return Value of the counting function.
    */
    virtual double count_fct(double r, double rc) const = 0;
    /**
    * Derivative of the counting function w.r.t. the distance.
    *
    * @param r distance between the two atoms 
    * @param rc summed up covalent radii of the two atoms
    * @return Derivative of the counting function.
    */
    virtual double dr_count_fct(double r, double rc) const = 0;
    /**
    * TCutoff function for large coordination numbers
    *
    * @param cn_max Maximum CN (not strictly obeyed).
    * @param cn On input coordination number, on output modified CN.
    * @param dcndr On input derivative of CN w.r.t. cartesian coordinates,
    * on output derivative of modified CN.
    * @param lgrad Flag for gradient calculation.
    */
    virtual int cut_coordination_number(const double cn_max, TVector<double> &cn, TMatrix<double> &dcndr, bool lgrad) = 0;
    // Constructor
    NCoordBase(double optional_kcn = 7.5, double optional_norm_exp = 1.0, double optional_cutoff = 25.0,
    double optional_f_directed = 1.0, double optional_cn_max = 8.0, const double* optional_rcov = rad);
    // Virtual destructor
    virtual ~NCoordBase() {
      cn.DelVec();
      dcndr.DelMat();
    }
};

// derived CN-class for erf()-based CN; includes default CN-cutoff
class NCoordErf : public NCoordBase {
  public:
    // erf() based counting function
    double count_fct(double, double) const override;
    // derivative of the erf() based counting function
    double dr_count_fct(double, double) const override;
    // Soft maximum/cutoff for coordination number
    int cut_coordination_number(const double, TVector<double>&, TMatrix<double>&, bool) override;
    // Constructor
    NCoordErf(double optional_kcn = 7.5, double optional_norm_exp = 1.0, double optional_cutoff = 25.0,
    double optional_f_directed = 1.0, double optional_cn_max = 8.0, const double* optional_rcov = rad)
    : NCoordBase(optional_kcn, optional_norm_exp, optional_cutoff, optional_f_directed, optional_cn_max, optional_rcov){}
    // Use default destructor; base class handles cleanup
    ~NCoordErf() override = default;
};

// derived CN-class for EEQ-BC local charges
class NCoordErfEN : public NCoordBase {
  public:
    // erf() based counting function
    double count_fct(double, double) const override;
    // derivative of the erf() based counting function
    double dr_count_fct(double, double) const override;
    // Soft maximum/cutoff for coordination number
    int cut_coordination_number(const double, TVector<double>&, TMatrix<double>&, bool) override;
    // coordination number scaling factor based on electronegativity difference
    double get_en_factor(int, int) const override;
    // Calculate the element pair-specific covalent radii for EEQ-BC qloc
    double get_rcov(int, int) const override;
    // Constructor
    NCoordErfEN(double optional_kcn = 7.5, double optional_norm_exp = 1.0, double optional_cutoff = 25.0,
    double optional_f_directed = -1.0, double optional_cn_max = 8.0, const double* optional_rcov = multicharge_param::eeqbc::eeqbc_cov_radii)
    : NCoordBase(optional_kcn, optional_norm_exp, optional_cutoff, optional_f_directed, optional_cn_max, optional_rcov){}
    // Use default destructor; base class handles cleanup
    ~NCoordErfEN() override = default;
};

// derived CN-class for the D4 model
class NCoordErfD4 : public NCoordBase {
  public:
    // erf() based counting function
    double count_fct(double, double) const override;
    // derivative of the erf() based counting function
    double dr_count_fct(double, double) const override;
    // coordination number scaling factor based on electronegativity difference
    double get_en_factor(int, int) const override;
    // Soft maximum/cutoff for coordination number
    int cut_coordination_number(const double, TVector<double>&, TMatrix<double>&, bool) override;
    // Constructor
    NCoordErfD4(double optional_kcn = 7.5, double optional_norm_exp = 1.0, double optional_cutoff = 25.0,
    double optional_f_directed = 1.0, double optional_cn_max = 8.0, const double* optional_rcov = rad)
    : NCoordBase(optional_kcn, optional_norm_exp, optional_cutoff, optional_f_directed, optional_cn_max, optional_rcov){}
    // Use default destructor; base class handles cleanup
    ~NCoordErfD4() override = default;
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

extern inline double log_cn_cut(double cn_max, double cn);

extern inline double dlog_cn_cut(double cn_max, double cn);

} // namespace dftd4
