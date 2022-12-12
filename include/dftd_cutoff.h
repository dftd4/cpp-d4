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

namespace dftd {

// Collection of real space cutoffs.
class TCutoff {
  public:
    double disp2;
    double disp3;
    double cn;
    double cn_eeq;

    explicit TCutoff();

    /**
     * @brief Set the cutoff for the two-body dispersion term.
     * 
     * @param val New value for two-body dispersion cutoff.
     */
    void set_disp2(double val);

    /**
     * @brief Set the cutoff for the three-body dispersion term.
     * 
     * @param val New value for three-body dispersion cutoff.
     */
    void set_disp3(double val);

    /**
     * @brief Set the cutoff for the DFT-D4 coordination number.
     * 
     * @param val New value for DFT-D4 coordination number cutoff.
     */
    void set_cn(double val);

    /**
     * @brief Set the cutoff for the EEQ coordination number.
     * 
     * @param val New value for EEQ coordination number cutoff.
     */
    void set_cn_eeq(double val);

    /**
     * @brief Set all cutoffs.
     * 
     * @param new_cutoff New value for all cutoffs (int/double).
     */
    void set_all(int new_cutoff);

    /**
     * @brief Set all cutoffs to the given value.
     * 
     * @param new_cutoff New value for all cutoffs (int/double).
     */
    void set_all(double new_cutoff);
};

}; // namespace dftd