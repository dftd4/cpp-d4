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
#ifndef TEST_GHOST_H
#define TEST_GHOST_H

#include <dftd_matrix.h>

using namespace dftd4;

static const double water_dimer_ref_cn[3]{
  +1.6104087812936088,
  +0.80554651093269469,
  +0.80486227036091429,
};

static const double water_dimer_ref_q[3]{
  -0.59265719057782107,
  +0.29739006625678255,
  +0.29526712432103869,
};

static const double water_dimer_ref_grad[3 * 6]{
  +3.6769324466783607E-05,
  +5.8683292759696172E-05,
  +0.0000000000000000E+00,
  +4.8209580157792990E-06,
  -4.4217699200268973E-05,
  +0.0000000000000000E+00,
  -4.1590282482562915E-05,
  -1.4465593559427221E-05,
  +0.0000000000000000E+00,
  // ghosts
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
  +0.0000000000000000E+00,
};

static const double water_ghost_ref_grad[3 * 6]{
  +0.0000000000000000E+00, // ghost
  +0.0000000000000000E+00, // ghost
  +0.0000000000000000E+00, // ghost
  +3.6769324466783607E-05, // O
  +5.8683292759696172E-05, // O
  +0.0000000000000000E+00, // O
  +0.0000000000000000E+00, // ghost
  +0.0000000000000000E+00, // ghost
  +0.0000000000000000E+00, // ghost
  +4.8209580157792990E-06, // H1
  -4.4217699200268973E-05, // H1
  +0.0000000000000000E+00, // H1
  -4.1590282482562915E-05, // H2
  -1.4465593559427221E-05, // H2
  +0.0000000000000000E+00, // H2
  +0.0000000000000000E+00, // ghost
  +0.0000000000000000E+00, // ghost
  +0.0000000000000000E+00, // ghost
};

extern int test_water(
  const int n,
  const char atoms[][4],
  const double coord[],
  const double ref_grad[],
  const TIVector &realIdx
);

extern int test_ghost(void);

#endif // TEST_GHOST_H
