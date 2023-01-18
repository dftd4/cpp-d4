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

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd4 {

// Default weighting factor for coordination number interpolation
static const double wf_default = 6.0;

// Default maximum charge scaling height for partial charge extrapolation
static const double ga_default = 3.0;

// Default charge scaling steepness for partial charge extrapolation
static const double gc_default = 2.0;

class TD4Model {
    public:
  double wf;
  double ga;
  double gc;

  explicit TD4Model(
    double wf_scale = wf_default,
    double ga_scale = ga_default,
    double gc_scale = gc_default
  );

  int weight_references(
    const TMolecule &mol,
    const TVector<double> &cn,
    const TVector<double> &q,
    TMatrix<double> &gwvec,
    TMatrix<double> &dgwdcn,
    TMatrix<double> &dgwdq,
    bool lgrad = false
  );

  int get_atomic_c6(
    const TMolecule &mol,
    const TMatrix<double> &gwvec,
    const TMatrix<double> &dgwdcn,
    const TMatrix<double> &dgwdq,
    TMatrix<double> &c6,
    TMatrix<double> &dc6dcn,
    TMatrix<double> &dc6dq,
    bool lgrad = false
  );

  int set_refalpha_eeq(const TMolecule &mol, TMatrix<double> &alpha);
};

extern inline double trapzd(const double a[23], const double b[23]);

extern inline double
  weight_cn(const double wf, const double cn, const double cnref);

extern inline double
  zeta(const double a, const double c, const double qref, const double qmod);
extern inline double
  dzeta(const double a, const double c, const double qref, const double qmod);

extern int get_max_ref(const TMolecule &mol, int &mref);

extern bool is_exceptional(double val);

} // namespace dftd4