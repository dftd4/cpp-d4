/* This file is part of cpp-d4.
 *
 * Copyright (C) 2019-2021 Sebastian Ehlert
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

// convert bohr (a.u.) to Ångström and back
static const double autoaa = 0.52917726;
static const double aatoau = 1.0 / autoaa;
// convert Hartree to eV and back
static const double autoev = 27.21138505;
static const double evtoau = 1.0 / autoev;
// convert Hartree to kcal/mol and back
static const double autokcal = 627.50947428;
static const double kcaltoau = 1.0 / autokcal;
// convert Hartree to kJ/mol and back
static const double autokj = 2625.49964038;
static const double kjtoau = 1.0 / autokj;
// convert Hartree to reciproce centimeters/wavenumbers and back
static const double autorcm = 219474.63067;
static const double autowav = autorcm;
static const double rcmtoau = 1.0 / autorcm;
static const double wavtoau = 1.0 / autowav;
// convert Hartree to nanometers and back
static const double autonm = 45.56335266;
static const double nmtoau = 1.0 / autonm;
// masses
// amu -> kg :: conversion from atomic mass units to kg
// me  -> kg :: electron mass (a.u.) in kg
// amu -> au :: conversion from a.u. to amu
static const double amutokg = 1.660539040e-27;
static const double kgtoamu = 1.0 / amutokg;
static const double metokg = 9.10938356e-31;
static const double kgtome = 1.0 / metokg;
static const double amutoau = amutokg * kgtome;
static const double autoamu = kgtoamu * metokg;
// femtosectons to atomic time units
static const double fstoau = 41.3413733365614;
static const double autofs = 1.0 / fstoau;
// Coulomb to atomic charge units (electrons)
static const double autoc = 1.6021766208e-19;
static const double ctoau = 1.0 / autoc;
