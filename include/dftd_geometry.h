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
  TMatrix<double> xyz;  // Cartesian Coordinates: (NAtoms x 3)-matrix
  TVector<int> at;   // atomic numbers

  TMolecule() { NAtoms = 0; }
  ~TMolecule() { FreeMemory(); }

  void GetMemory(int NumAt_) {
    FreeMemory();
    NAtoms = NumAt_;
    xyz.New(NAtoms, 3);
    at.New(NAtoms);
  }

  void FreeMemory(void) {
    xyz.Delete();
    at.Delete();
  }

  /* set structure by array of coords (double), atom types (int) and number of atoms */
  /* Geometry has to be in bohrs */
  void SetGeometry(const double* coord, const int* attype, int number)
  {
      GetMemory(number);
      for(int i = 0; i < NAtoms; ++i)
      {
          at(i) = attype[i];
          xyz(i, 0) = coord[3*i + 0];
          xyz(i, 1) = coord[3*i + 1];
          xyz(i, 2) = coord[3*i + 2];
      }
  }

  /* update the geometry */
  /* Geometry has to be in bohrs */
  void UpdateGeometry(const double* coord)
  {
      for(int i = 0; i < NAtoms; ++i)
      {
          xyz(i, 0) = coord[3*i + 0];
          xyz(i, 1) = coord[3*i + 1];
          xyz(i, 2) = coord[3*i + 2];
      }
  }

  /* update the position of a single atom (index, coord and element) */
  /* Geometry has to be in bohrs */
  void UpdateAtom(int index, double x, double y, double z, int element = -1)
  {
      //std::cout << x << " " << y << " " << z << std::endl;
      if(element != -1)
        at(index) = element;
      xyz(index, 0) = x;
      xyz(index, 1) = y;
      xyz(index, 2) = z;
  }

  void Print() const
  {
      for(int i = 0; i < NAtoms; ++i)
      {
          std::cout << at(i) << " " << xyz(i, 0) << " " << xyz(i, 1) << " " << xyz(i, 2) << std::endl;

      }
  }
};

}  // namespace dftd4
