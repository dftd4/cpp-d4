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

#include "dftd_readxyz.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

#include "dftd_econv.h"
#include "dftd_geometry.h"

void read_xyzfile(const std::string& name, dftd::TMolecule& mol) {
  std::ifstream geo;
  std::stringstream buffer;
  std::string line;
  int n{0};
  int at{0};
  char sym[3]{"  "};
  double x{0.0}, y{0.0}, z{0.0};

  geo.open(name);
  if (!geo) {
    printf("Error: Cannot open file.");
    exit(EXIT_FAILURE);
  } 

  buffer << (std::getline(geo, line), line);
  buffer >> n;

  mol.GetMemory(n);

  std::getline(geo, line);  // skip comment line

  for (int i = 0; i != n; i++) {
    geo >> sym >> x >> y >> z;
    if (geo.fail()) {
      printf("Error: Input file could not be read.");
      exit(EXIT_FAILURE);
    }
    mol.xyz(i, 0) = x * aatoau;
    mol.xyz(i, 1) = y * aatoau;
    mol.xyz(i, 2) = z * aatoau;
    at = element(sym);
    mol.at(i) = at;
  }

  geo.close();
}

int element(const std::string& sym) {
  char elem[3]{"  "};
  char pse[118][3]{
      "h ", "he", "li", "be", "b ", "c ", "n ", "o ", "f ", "ne", "na", "mg",
      "al", "si", "p ", "s ", "cl", "ar", "k ", "ca", "sc", "ti", "v ", "cr",
      "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr",
      "rb", "sr", "y ", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd",
      "in", "sn", "sb", "te", "i ", "xe", "cs", "ba", "la", "ce", "pr", "nd",
      "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf",
      "ta", "w ", "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po",
      "at", "rn", "fr", "ra", "ac", "th", "pa", "u ", "np", "pu", "am", "cm",
      "bk", "cf", "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs",
      "mt", "ds", "rg", "cn", "nh", "fl", "mc", "lv", "ts", "og",
  };

  std::transform(sym.begin(), sym.end(), elem, ::tolower);

  for (int i = 0; i != 118; i++) {
    std::string a = elem;
    std::string b = pse[i];
    int stat = a.compare(b);
    if (!stat) return i + 1;
  }

  return EXIT_SUCCESS;
}
