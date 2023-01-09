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
#include "dftd_cutoff.h"

namespace dftd4 {

TCutoff::TCutoff(
  double cut_disp2 /* = disp2_default*/,
  double cut_disp3 /* = disp3_default*/,
  double cut_cn /* = cn_default*/,
  double cut_cn_eeq /* = cn_eeq_default*/
) {
  disp2 = cut_disp2;
  disp3 = cut_disp3;
  cn = cut_cn;
  cn_eeq = cut_cn_eeq;
}

} // namespace dftd4
