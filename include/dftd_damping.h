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

/**
 * Collection of damping parameters for D4-BJ-EEQ-ATM and D4-BJ-EEQ-MBD.
 * @see https://github.com/dftd4/dftd4/blob/main/assets/parameters.toml
 */
#pragma once

#include <string>

#include "dftd_dispersion.h"

namespace dftd4 {

/**
 * @brief Collect the D4 parameters.
 *
 * @param func Name of the functional.
 * @param par Parameter class containing the parameters for selected functional.
 * @param latm Switch for D4-ATM (true) or D4-MBD (false) parameters.
 * @returns Exit status.
 */
extern int
  d4par(const std::string func, dftd4::dparam &par, bool latm = true);

} // namespace dftd4
