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

#include <string>

#include "dftd_geometry.h"

static const char pse[118][3]{
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

extern void read_xyzfile(const std::string &, dftd4::TMolecule &);

extern int element(const std::string &);
