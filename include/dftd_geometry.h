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

#pragma once

#include "dftd_matrix.h"

namespace dftd4 {
// Input of the molecular geometry
class TMolecule {
    public:
  int NAtoms;
  TMatrix<double> CC; // Cartesian Coordinates: (NAtoms x 3)-matrix
  TVector<int> ATNO;     // atomic numbers

  TMolecule() { NAtoms = 0; }
  ~TMolecule() { FreeMemory(); }

  void GetMemory(int NumAt_) {
    FreeMemory();
    NAtoms = NumAt_;
    CC.New(NAtoms, 3);
    ATNO.New(NAtoms);
  }

  void FreeMemory(void) {
    CC.Delete();
    ATNO.Delete();
  }
};

} // namespace dftd4
