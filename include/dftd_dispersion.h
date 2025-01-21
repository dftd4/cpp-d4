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
 * D4(EEQ)-ATM implementation
 */
#pragma once

#include "dftd_cutoff.h"
#include "dftd_geometry.h"
#include "dftd_model.h"

namespace dftd4 {

/**
 * @brief Container for DFT-D4 damping parameters.
 *
 * @param s6 C6 damping parameter. Usually 1.
 * @param s8 C8 damping parameter.
 * @param s10 C10 damping parameter. Always 0 currently.
 * @param s9 ATM (three-body) damping parameter. Usually 1.
 * @param a1 Becke-Johnson damping parameter.
 * @param a2 Becke-Johnson damping parameter.
 * @param alp ATM (three-body) damping parameter. Usually 16.
 */
class dparam {
  public:
    double s6;
    double s8;
    double s10;
    double s9;
    double a1;
    double a2;
    int alp;
};

/**
 * @brief Wrapper to handle the evaluation of dispersion energy and derivatives.
 *
 * @param mol Molecular geometry.
 * @param charge Molecular charge.
 * @param par DFT-D4 parameters.
 * @param d4 Base D4 dispersion model.
 * @param cutoff Real-space cutoffs for CN and dispersion.
 * @param energy Dispersion energy (inout).
 * @param GRAD Dispersion gradient (inout).
 * @return Exit status.
 */
extern int get_dispersion(
  const TMolecule &mol,
  int charge,
  const TD4Model &d4,
  const dparam &par,
  TCutoff cutoff,
  double &energy,
  double *GRAD
);

/**
 * @brief Wrapper to handle the evaluation of dispersion energy and derivatives.
 *
 * This function calculates the atom-wise dispersion energy and gradients 
 * for the given molecular geometry, considering only the atoms specified 
 * in `realIdx`.
 *
 * @param mol Molecular geometry.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param charge Molecular charge.
 * @param par DFT-D4 parameters.
 * @param d4 Base D4 dispersion model.
 * @param cutoff Real-space cutoffs for CN and dispersion.
 * @param energies atom-wise dispersion energies (inout).
 * @param GRAD Dispersion gradient (inout).
 * @return Exit status.
 */
extern int get_dispersion(
  const TMolecule &mol,
  const TIVector &realIdx,
  int charge,
  const TD4Model &d4,
  const dparam &par,
  TCutoff cutoff,
  TRVector &energies,
  double *GRAD
);

/**
 * @brief Wrapper to handle the evaluation of dispersion energy and derivatives.
 *
 * This function calculates the dispersion energy and gradients for the given
 * molecular geometry, considering only the atoms specified in `realIdx`.
 *
 * @param mol Molecular geometry.
 * @param realIdx List for real atoms excluding ghost/non atoms.
 * @param charge Molecular charge.
 * @param par DFT-D4 parameters.
 * @param d4 Base D4 dispersion model.
 * @param cutoff Real-space cutoffs for CN and dispersion.
 * @param energy Dispersion energy (inout).
 * @param GRAD Dispersion gradient (inout).
 * @return Exit status.
 */
extern int get_dispersion(
  const TMolecule &mol,
  const TIVector &realIdx,
  int charge,
  const TD4Model &d4,
  const dparam &par,
  TCutoff cutoff,
  double &energy,
  double *GRAD
);

} // namespace dftd4
