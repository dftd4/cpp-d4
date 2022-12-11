/* This file is part of cpp-d4.
 *
 * Copyright (C) 2019-2021 Sebastian Ehlert
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

namespace dftd {

/**
 * Calculate all distance pairs and store in matrix
 * 
 * @param mol Molecule object.
 * @param dist Distance matrix (inout).
 * @return Exit status.
*/
extern int calc_distances(TMolecule& mol, TMatrix<double>& dist);

/**
 * Calculate error function coordination number.
 * 
 * @param mol Molecule object.
 * @param dist Distance matrix.
 * @param cn Vector of coordination numbers.
 * @param thr Real-space cutoff (default: 1600.0).
 * @return Exit status.
*/
extern int ncoord_erf(TMolecule& mol, TMatrix<double>& dist,
                      TVector<double>& cn, double thr = 1600.0);

/**
 * Calculate covalent coordination number for DFT-D4.
 * 
 * @param mol Molecule object.
 * @param dist Distance matrix.
 * @param cn Vector of coordination numbers.
 * @param thr Real-space cutoff (default: 1600.0).
 * @return Exit status.
*/
extern int ncoord_d4(TMolecule& mol, TMatrix<double>& dist,
                     TVector<double>& cn, double thr = 1600.0);

/**
 * Calculate error function coordination number and derivative
 * w.r.t. nuclear coordinates.
 * 
 * @param mol Molecule object.
 * @param dist Distance matrix.
 * @param cn Vector of coordination numbers.
 * @param dcndr Derivative of coordination number.
 * @param thr Real-space cutoff (default: 1600.0).
 * @return Exit status.
*/
extern int dncoord_erf(TMolecule& mol, TMatrix<double>& dist,
                       TVector<double>& cn, TMatrix<double>& dcndr,
                       double thr = 1600.0);

/**
 * Calculate covalent coordination number and derivative
 * w.r.t. nuclear coordinates
 * 
 * @param mol Molecule object.
 * @param dist Distance matrix.
 * @param cn Vector of coordination numbers.
 * @param dcndr Derivative of coordination number.
 * @param thr Real-space cutoff (default: 1600.0).
 * @return Exit status.
*/
extern int dncoord_d4(TMolecule& mol, TMatrix<double>& dist,
                      TVector<double>& cn, TMatrix<double>& dcndr,
                      double thr = 1600.0);

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
 * @param rr Cutoff radius.
 * @return Derivative of the counting function.
 */
extern double derf_count(double k, double rr);

};  // namespace dftd
