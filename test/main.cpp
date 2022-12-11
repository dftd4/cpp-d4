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
#include <cstdlib>
#include <iostream>
#include <map>

#include "test_disp.h"
#include "test_ncoord.h"


enum test {
  invalid,
  disp,
  ncoord,
};

test get_tests(std::string name) {
  static const std::map<std::string, test> testStrings {
      {"disp", disp},
      {"ncoord", ncoord},
  };

  std::string test = name;
  transform(name.begin(), name.end(), test.begin(), ::tolower);
  auto iter = testStrings.find(test);
  if (iter != testStrings.end()) {
    return iter->second;
  }
  return invalid;
};


int main(int argc, char *argv[]) {
  switch (get_tests(argv[1])) {
    default:
      return EXIT_FAILURE;
    case ncoord:
      return test_ncoord(); 
    case disp:
      return test_disp(); 
  }
}
