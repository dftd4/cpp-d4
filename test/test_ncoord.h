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
#ifndef TEST_NCOORD_H
#define TEST_NCOORD_H

static const double mb16_43_01_ref_cn[16] {
  +3.0734967711040162E+00, +9.3146160511610288E-01, +1.4370943937583889E+00, 
  +1.3330943158196045E+00, +7.2074352703033706E-01, +8.5965900477098189E-01, 
  +1.3578215817792136E+00, +1.5394000699602508E+00, +3.1940036819525854E+00, 
  +8.1216211163134222E-01, +8.5953344378485375E-01, +1.5334710815558710E+00, 
  +4.2331498952572053E+00, +3.0304850456739563E+00, +3.4522931948830609E+00, 
  +4.2847828965226427E+00
};

static const double rost61_m1_ref_cn[22] {
  +3.7743890839860823E+00, +3.7533525980157583E+00, +3.7272213911500702E+00, 
  +9.2591719422049701E-01, +3.7403891707352228E+00, +3.7389943427981649E+00, 
  +9.2590197299758326E-01, +9.2599894955767703E-01, +9.2611826399971675E-01, 
  +9.2576484754843247E-01, +3.7170298017526591E+00, +3.7298933986021074E+00, 
  +3.7649883287219090E+00, +9.2586986693780782E-01, +3.7529905532199597E+00, 
  +3.7621775881176829E+00, +9.2587667173548660E-01, +9.2604216237640058E-01, 
  +9.2598322975544489E-01, +9.2588963692121407E-01, +8.1580076383956079E+00, 
  +8.4210289216027034E-01
};

extern int test_cn(
  int n,
  const char atoms[][3],
  const double coord[],
  const double ref_cn[]
);

extern int test_ncoord();

#endif // TEST_NCOORD_H
