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
#ifndef TEST_GRAD_H
#define TEST_GRAD_H

#include "dftd_dispersion.h"
#include "molecules.h"

using namespace dftd4;

extern int test_numgrad(TMolecule &mol, const int charge, const dparam &par);
extern int is_trans_invar(const TMolecule &mol, double gradient[]);

extern int test_bp86d4atm_water();
extern int test_pbed4_mb01();
extern int test_tpss0d4mbd_rost61m1();

extern int test_grad();

#endif // TEST_GRAD_H
