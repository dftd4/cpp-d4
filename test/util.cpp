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
#include <algorithm>
#include <cmath>
#include <map>
#include <string>

#include <dftd_geometry.h>
#include <dftd_matrix.h>

#include "util.h"

using namespace dftd4;

int get_molecule(
  int n,
  const char atoms[][3],
  const double coord[],
  TMolecule &mol
) {
  ;
  mol.GetMemory(n);
  for (int i = 0; i != n; i++) {
    mol.CC(i, 0) = coord[3 * i];
    mol.CC(i, 1) = coord[3 * i + 1];
    mol.CC(i, 2) = coord[3 * i + 2];
    mol.ATNO(i) = element(atoms[i]);
  }

  return EXIT_SUCCESS;
}

bool check(
  double actual,
  double expected,
  double epsilon /*= 1e-12*/,
  bool rel /*= false*/
) {
  double diff;

  if (rel) {
    diff = std::fabs(actual - expected) / expected;
  } else {
    diff = std::fabs(actual - expected);
  }

  if (diff > epsilon) { return EXIT_FAILURE; }
  return EXIT_SUCCESS;
};

bool check(
  float actual,
  float expected,
  float epsilon /*= 1e-6*/,
  bool rel /*= false*/
) {
  float diff;

  if (rel) {
    diff = std::fabs(actual - expected) / expected;
  } else {
    diff = std::fabs(actual - expected);
  }

  if (diff > epsilon) { return EXIT_FAILURE; }
  return EXIT_SUCCESS;
};

void print_fail(const char specifier[32], double obtained, double expected) {
  printf("Failed for: '%s'\n", specifier);
  printf("  - Expected   %.16f\n", expected);
  printf("  - Obtained   %.16f\n", obtained);
  printf("  - Diff (abs)  %.16f\n", fabs(expected - obtained));
}

int element(const std::string &sym) {
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
