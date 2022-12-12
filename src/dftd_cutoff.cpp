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

namespace dftd {

// Real space cutoff for CN within D4
static const double cn_default = 30.0;

// Real space cutoff for CN within EEQ
static const double cn_eeq_default = 25.0;

// Two-body interaction cutoff
static const double disp2_default = 60.0;

// Three-body interaction cutoff
static const double disp3_default = 40.0;


TCutoff::TCutoff() {
  disp2 = disp2_default;
  disp3 = disp3_default;
  cn = cn_default;
  cn_eeq = cn_eeq_default;
};

void TCutoff::set_all(int new_cutoff) {
  disp2 = (double)new_cutoff;
  disp3 = (double)new_cutoff;
  cn = (double)new_cutoff;
  cn_eeq = (double)new_cutoff;
};
void TCutoff::set_all(double new_cutoff) {
  disp2 = new_cutoff;
  disp3 = new_cutoff;
  cn = new_cutoff;
  cn_eeq = new_cutoff;
};

void TCutoff::set_disp2(double val) {
  disp2 = val;
};

void TCutoff::set_disp3(double val) {
  disp3 = val;
};

void TCutoff::set_cn(double val) {
  cn = val;
};

void TCutoff::set_cn_eeq(double val) {
  cn_eeq = val;
};

}; // namespace dftd