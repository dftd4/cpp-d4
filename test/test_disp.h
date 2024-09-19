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
#ifndef TEST_DISP_H
#define TEST_DISP_H

static const double water_ref_energy{-2.3162148796416E-04};

static const double mb16_43_01_ref_energy{-2.5882588037023E-02};

static const double rost61_m1_ref_energy{-3.4287391104745E-02};

static const double amf3_ref_energy{-2.9119716615478E-03};

static const double actinides_ref_energy{-2.836190476292469E-01};

extern int test_energy(
  int n,
  const char atoms[][4],
  const double coord[],
  int charge,
  double ref
);

extern int test_disp();

#endif // TEST_DISP_H
