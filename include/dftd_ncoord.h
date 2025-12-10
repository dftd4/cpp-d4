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
#include <cmath>

namespace dftd4 {

class NCoordBase {
  public:
    static const double
      rad[];    // cov. radii default: D3 covalent radii from dftd_ncoord.cpp
    double kcn; // Steepness of counting function
    double
      norm_exp; // exponent of the normalizing-factor in the counting function
    double cutoff;      // Coordination number cutoff distance
    double f_directed;  // directed factor for EN scaled coordination number
    double cn_max;      // (Soft-)maximum value for the coordination number
    const double *rcov; // covalent radii for CN
    // Get the coordination number
    /**
     * Wrapper for coordination number calculation.
     *
     * @param mol Molecule object.
     * @param dist Distance matrix.
     * @param cn Coordination number.
     * @param dcndr Derivative of the coordination number.
     * @param lgrad Flag for gradient computation.
     * @return Exit status.
     */
    int get_ncoord( // without ghost atom indices
      const TMolecule &mol,
      const TMatrix<double> &dist,
      TVector<double> &cn,  // coordination number
      TMatrix<double> &dcndr,  // derivative of the coordination number
      bool lgrad);
    /**
     * Wrapper for coordination number calculation.
     * Allocates the cn vector unconditionally, and the dcndr matrix only if
     * lgrad is true.
     *
     * @param mol Molecule object.
     * @param realIdx List for real atoms excluding ghost/non atoms.
     * @param dist Distance matrix.
     * @param cn Coordination number.
     * @param dcndr Derivative of the coordination number.
     * @param lgrad Flag for gradient computation.
     * @return Exit status.
     */
    int get_ncoord(  // with ghost atoms
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TVector<double> &cn,  // coordination number
      TMatrix<double> &dcndr,  // derivative of the coordination number
      bool lgrad);
    /**
     * Calculate the coordination number.
     *
     * @param mol Molecule object.
     * @param realIdx List for real atoms excluding ghost/non atoms.
     * @param dist Distance matrix.
     * @param cn Coordination number.
     * @param dcndr Derivative of the coordination number.
     * @return Exit status.
     */
    int ncoord_base(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TVector<double> &cn
    );
    /**
     * Calculate error function coordination number and derivative
     * w.r.t. nuclear coordinates.
     *
     * @param mol Molecule object.
     * @param realIdx List for real atoms excluding ghost/non atoms.
     * @param dist Distance matrix.
     * @param cn Coordination number.
     * @param dcndr Derivative of the coordination number.
     * @return Exit status.
     */
    int dr_ncoord_base(
      const TMolecule &mol,
      const TIVector &realIdx,
      const TMatrix<double> &dist,
      TVector<double> &cn,
      TMatrix<double> &dcndr
    );
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
    virtual int cut_coordination_number(
      const double cn_max,
      TVector<double> &cn,
      TMatrix<double> &dcndr,
      bool lgrad
    ) = 0;
    // Constructor
    NCoordBase(
      double optional_kcn = 7.5,
      double optional_norm_exp = 1.0,
      double optional_cutoff = 25.0,
      double optional_f_directed = 1.0,
      double optional_cn_max = 8.0,
      const double *optional_rcov = rad
    );
    // Virtual destructor
    virtual ~NCoordBase() {}
};

inline double NCoordBase::get_en_factor(int i, int j) const { return 1.0; }

inline double NCoordBase::get_rcov(int i, int j) const {
  return rcov[i] + rcov[j];
}

// derived CN-class for erf()-based CN; includes default CN-cutoff
class NCoordErf : public NCoordBase {
  public:
    // erf() based counting function
    double count_fct(double, double) const override;
    // derivative of the erf() based counting function
    double dr_count_fct(double, double) const override;
    // Soft maximum/cutoff for coordination number
    int cut_coordination_number(
      const double,
      TVector<double> &,
      TMatrix<double> &,
      bool
    ) override;
    // Constructor
    NCoordErf(
      double optional_kcn = 7.5,
      double optional_norm_exp = 1.0,
      double optional_cutoff = 25.0,
      double optional_f_directed = 1.0,
      double optional_cn_max = 8.0,
      const double *optional_rcov = rad
    )
      : NCoordBase(
          optional_kcn,
          optional_norm_exp,
          optional_cutoff,
          optional_f_directed,
          optional_cn_max,
          optional_rcov
        ) {}
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
    int cut_coordination_number(
      const double,
      TVector<double> &,
      TMatrix<double> &,
      bool
    ) override;
    // coordination number scaling factor based on electronegativity difference
    double get_en_factor(int, int) const override;
    // Calculate the element pair-specific covalent radii for EEQ-BC qloc
    double get_rcov(int, int) const override;
    // Constructor
    NCoordErfEN(
      double optional_kcn = 7.5,
      double optional_norm_exp = 1.0,
      double optional_cutoff = 25.0,
      double optional_f_directed = -1.0,
      double optional_cn_max = 8.0,
      const double *optional_rcov = multicharge_param::eeqbc::eeqbc_cov_radii
    )
      : NCoordBase(
          optional_kcn,
          optional_norm_exp,
          optional_cutoff,
          optional_f_directed,
          optional_cn_max,
          optional_rcov
        ) {}
    // Use default destructor; base class handles cleanup
    ~NCoordErfEN() override = default;
};

inline double NCoordErfEN::get_en_factor(int i, int j) const {
  return multicharge_param::eeqbc::eeqbc_en[j] -
         multicharge_param::eeqbc::eeqbc_en[i];
}

inline double NCoordErfEN::get_rcov(int i, int j) const {
  return rcov[i] + rcov[j];
}

// derived CN-class for the D4 model
class NCoordErfD4 : public NCoordBase {
  private:
    static constexpr double k4 = 4.10451;
    static constexpr double k5 = 19.08857;
    static constexpr double k6 = 2 * 11.28174 * 11.28174;

  public:
    // erf() based counting function
    double count_fct(double, double) const override;
    // derivative of the erf() based counting function
    double dr_count_fct(double, double) const override;
    // coordination number scaling factor based on electronegativity difference
    double get_en_factor(int, int) const override;
    // Soft maximum/cutoff for coordination number
    int cut_coordination_number(
      const double,
      TVector<double> &,
      TMatrix<double> &,
      bool
    ) override;
    // Constructor
    NCoordErfD4(
      double optional_kcn = 7.5,
      double optional_norm_exp = 1.0,
      double optional_cutoff = 25.0,
      double optional_f_directed = 1.0,
      double optional_cn_max = 8.0,
      const double *optional_rcov = rad
    )
      : NCoordBase(
          optional_kcn,
          optional_norm_exp,
          optional_cutoff,
          optional_f_directed,
          optional_cn_max,
          optional_rcov
        ) {}
    // Use default destructor; base class handles cleanup
    ~NCoordErfD4() override = default;
};

// pauling EN's
static const double en[119]{
  0.0,  // dummy
  2.20, // H
  3.00, // He
  0.98, // Li (2nd)
  1.57, // Be
  2.04, // B
  2.55, // C
  3.04, // N
  3.44, // O
  3.98, // F
  4.50, // Ne
  0.93, // Na (3rd)
  1.31, // Mg
  1.61, // Al
  1.90, // Si
  2.19, // P
  2.58, // S
  3.16, // Cl
  3.50, // Ar
  0.82, // K  (4th)
  1.00, // Ca
  1.36, // Sc
  1.54, // Ti
  1.63, // V
  1.66, // Cr
  1.55, // Mn
  1.83, // Fe
  1.88, // Co
  1.91, // Ni
  1.90, // Cu
  1.65, // Zn
  1.81, // Ga
  2.01, // Ge
  2.18, // As
  2.55, // Se
  2.96, // Br
  3.00, // Kr
  0.82, // Rb (5th)
  0.95, // Sr
  1.22, // Y
  1.33, // Zr
  1.60, // Nb
  2.16, // Mo
  1.90, // Tc
  2.20, // Ru
  2.28, // Rh
  2.20, // Pd
  1.93, // Ag
  1.69, // Cd
  1.78, // In
  1.96, // Sn
  2.05, // Sb
  2.10, // Te
  2.66, // I
  2.60, // Xe
  0.79, // Cs (6th)
  0.89, // Ba
  1.10, // La
  1.12, // Ce
  1.13, // Pr
  1.14, // Nd
  1.15, // Pm
  1.17, // Sm
  1.18, // Eu
  1.20, // Gd
  1.21, // Tb
  1.22, // Dy
  1.23, // Ho
  1.24, // Er
  1.25, // Tm
  1.26, // Yb
  1.27, // Lu
  1.30, // Hf
  1.50, // Ta
  2.36, // W
  1.90, // Re
  2.20, // Os
  2.20, // Ir
  2.28, // Pt
  2.54, // Au
  2.00, // Hg
  1.62, // Tl
  2.33, // Pb
  2.02, // Bi
  2.00, // Po
  2.20, // At
  2.20, // Rn
  0.79, // Fr (7th)
  0.90, // Ra
  1.10, // Ac
  1.30, // Th
  1.50, // Pa
  1.38, // U
  1.36, // Np
  1.28, // Pu
  1.30, // Am
  1.30, // Cm
  1.30, // Bk
  1.30, // Cf
  1.30, // Es
  1.30, // Fm
  1.30, // Md
  1.30, // No
  1.30, // Lr
  1.50, // Rf (only dummies from here)
  1.50, // Db
  1.50, // Sg
  1.50, // Bh
  1.50, // Hs
  1.50, // Mt
  1.50, // Ds
  1.50, // Rg
  1.50, // Cn
  1.50, // Nh
  1.50, // Fl
  1.50, // Lv
  1.50, // Mc
  1.50, // Ts
  1.50, // Og
};

inline double NCoordErfD4::get_en_factor(int i, int j) const {
  return k4 * exp(-pow((fabs(en[i] - en[j]) + k5), 2) / k6);
}

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
