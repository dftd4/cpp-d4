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
#include "test_disp2.h"
#include "test_ghost.h"
#include "test_grad.h"
#include "test_multi.h"
#include "test_ncoord.h"
#include "test_param.h"

enum test {
  invalid,
  disp2,
  disp,
  ghost,
  grad,
  ncoord,
  param,
  multi,
};

test get_tests(std::string name) {
  static const std::map<std::string, test> testStrings{
    {"disp", disp},
    {"disp2", disp2},
    {"ghost", ghost},
    {"grad", grad},
    {"ncoord", ncoord},
    {"param", param},
    {"multi", multi}
  };

  std::string test = name;
  transform(name.begin(), name.end(), test.begin(), ::tolower);
  auto iter = testStrings.find(test);
  if (iter != testStrings.end()) { return iter->second; }
  return invalid;
};

int main(int argc, char *argv[]) {
  switch (get_tests(argv[1])) {
  default: return EXIT_FAILURE;
  case disp2: return test_disp2();
  case disp: return test_disp();
  case ghost: return test_ghost();
  case grad: return test_grad();
  case ncoord: return test_ncoord();
  case param: return test_param();
  case multi: return test_multi();
  }
}
