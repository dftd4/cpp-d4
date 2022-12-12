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
 * Definition of the coordination numbers used in DFT-D4 and the
 * electronegativity equilibration (EEQ) model
 *
 * This module works on a distance matrix to avoid recalculating
 * the distances every time.
 */
#include "dftd_ncoord.h"

#include <cmath>
#include <iostream>

#include "dftd_geometry.h"
#include "dftd_matrix.h"

namespace dftd {

// convert bohr (a.u.) to Ångström and back
static const double autoaa = 0.52917721090449243;
static const double aatoau = 1.0 / autoaa;

/**
 * Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
 * 188-197), values for metals decreased by 10 %.
 * 
 * These values are actually never used in the code.
 * Only the scaled values below are used (`rad`).
*/
static const double covalent_rad_d3[119]{ 0.0,
    0.32, 0.46,                                                  //  H,He
    1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,              //  Li-Ne
    1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,              //  Na-Ar
    1.76, 1.54,                                                  //  K,Ca
    1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,  //  Sc-Zn
    1.12, 1.09, 1.15, 1.10, 1.14, 1.17,                          //  Ga-Kr
    1.89, 1.67,                                                  //  Rb,Sr
    1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23,  //  Y-Cd
    1.28, 1.26, 1.26, 1.23, 1.32, 1.31,                          //  In-Xe
    2.09, 1.76,                                                  //  Cs,Ba
    1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51,                    //  La-Eu
    1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,                    //  Gd-Yb
    1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,  //  Lu-Hg
    1.30, 1.30, 1.36, 1.31, 1.38, 1.42,                          //  Tl-Rn
    2.01, 1.81,                                                  //  Fr,Ra
    1.67, 1.58, 1.52, 1.53, 1.54, 1.55, 1.49,                    //  Ac-Am
    1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58,                    //  Cm-No
    1.45, 1.41, 1.34, 1.29, 1.27, 1.21, 1.16, 1.15, 1.09, 1.22,  //  Lr-Cn
    1.22, 1.29, 1.46, 1.58, 1.48, 1.41                           //  Nh-Og
};

/**
 * D3 covalent radii used to construct the coordination number.
 * rad = covalent_rad_d3 * 4.0/3.0 * aatoau
 */
static const double rad[119]{ 0.0,
    0.8062831465047213,     1.1590320231005369,     3.0235617993927044,
    2.3684567428576182,     1.9401188212769855,     1.8897261246204402,
    1.7889407313073502,     1.5873699446811698,     1.6125662930094427,
    1.6881553379942602,     3.5274887659581555,     3.1495435410340673,
    2.8471873610947966,     2.6204202261403440,     2.7715983161099791,
    2.5700275294837986,     2.4944384844989811,     2.4188494395141635,
    4.4345573057759662,     3.8802376425539711,     3.3511143276602473,
    3.0739544960492493,     3.0487581477209771,     2.7715983161099791,
    2.6960092711251620,     2.6204202261403440,     2.5196348328272538,
    2.4944384844989811,     2.5448311811555264,     2.7464019677817069,
    2.8219910127665244,     2.7464019677817069,     2.8975800577513415,
    2.7715983161099791,     2.8723837094230693,     2.9479727544078864,
    4.7621098340435095,     4.2077901708215135,     3.7038632042560633,
    3.5022924176298824,     3.3259179793319751,     3.1243471927057946,
    2.8975800577513415,     2.8471873610947966,     2.8471873610947966,
    2.7212056194534342,     2.8975800577513415,     3.0991508443775224,
    3.2251325860188853,     3.1747398893623395,     3.1747398893623395,
    3.0991508443775224,     3.3259179793319751,     3.3007216310037024,
    5.2660368006089602,     4.4345573057759662,     4.0818084291801515,
    3.7038632042560633,     3.9810230358670613,     3.9558266875387886,
    3.9306303392105164,     3.9054339908822433,     3.8046485975691535,
    3.8298449458974257,     3.8046485975691535,     3.7794522492408804,
    3.7542559009126082,     3.7542559009126082,     3.7290595525843355,
    3.8550412942256984,     3.6786668559277906,     3.4518997209733380,
    3.3007216310037024,     3.0991508443775224,     2.9731691027361595,
    2.9227764060796142,     2.7967946644382522,     2.8219910127665244,
    2.8471873610947966,     3.3259179793319751,     3.2755252826754302,
    3.2755252826754302,     3.4267033726450653,     3.3007216310037024,
    3.4770960693016097,     3.5778814626147004,     5.0644660139827797,
    4.5605390474173291,     4.2077901708215135,     3.9810230358670613,
    3.8298449458974257,     3.8550412942256984,     3.8802376425539711,
    3.9054339908822433,     3.7542559009126082,     3.7542559009126082,
    3.8046485975691535,     3.8046485975691535,     3.7290595525843355,
    3.7794522492408804,     3.9306303392105164,     3.9810230358670613,
    3.6534705075995175,     3.5526851142864277,     3.3763106759885204,
    3.2503289343471575,     3.1999362376906122,     3.0487581477209771,
    2.9227764060796142,     2.8975800577513415,     2.7464019677817069,
    3.0739544960492493,     3.4267033726450653,     3.6030778109429726,
    3.6786668559277906,     3.9810230358670613,     3.7290595525843355,
    3.9558266875387886 
};

// pauling EN's
static const double en[119]{
    0.0, 2.20, 3.00,                                             //  H,He
    0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 4.50,              //  Li-Ne
    0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 3.50,              //  Na-Ar
    0.82, 1.00,                                                  //  K,Ca
    1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65,  //  Sc-Zn
    1.81, 2.01, 2.18, 2.55, 2.96, 3.00,                          //  Ga-Kr
    0.82, 0.95,                                                  //  Rb,Sr
    1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69,  //  Y-Cd
    1.78, 1.96, 2.05, 2.10, 2.66, 2.60,                          //  In-Xe
    0.79, 0.89,                                                  //  Cs,Ba
    1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18,                    //  La-Eu
    1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26,                    //  Gd-Yb
    1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00,  //  Lu-Hg
    1.62, 2.33, 2.02, 2.00, 2.20, 2.20,                          //  Tl-Rn
    // only dummies below
    1.50, 1.50,                                                  //  Fr,Ra
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                    //  Ac-Am
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                    //  Cm-No
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,  //  Rf-Cn
    1.50, 1.50, 1.50, 1.50, 1.50, 1.50                           //  Nh-Og
};

static const double kn = 7.5;
static const double k4 = 4.10451;
static const double k5 = 19.08857;
static const double k6 = 2 * pow(11.28174, 2);
static const double hlfosqrtpi = 1.0 / 1.77245385091;
static const double cn_max = 8.0;


int calc_distances(const TMolecule& mol, TMatrix<double>& dist) {
  double rx = 0.0, ry = 0.0, rz = 0.0, tmp = 0.0;
  for (int i = 0; i != mol.NAtoms; i++) {
    dist(i, i) = 0.0;
    for (int j = 0; j != i; j++) {
      rx = mol.xyz(i, 0) - mol.xyz(j, 0);
      ry = mol.xyz(i, 1) - mol.xyz(j, 1);
      rz = mol.xyz(i, 2) - mol.xyz(j, 2);
      tmp = sqrt(rx * rx + ry * ry + rz * rz);
      dist(i, j) = tmp;
      dist(j, i) = tmp;
    }
  }

  return EXIT_SUCCESS;
}


int get_ncoord_erf(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  TVector<double>& cn,
  TMatrix<double>& dcndr,
  bool lgrad,
  double thr
) {
  int info;

  if (lgrad) {
    info = dncoord_erf(mol, dist, cn, dcndr, thr);
  } else {
    info = ncoord_erf(mol, dist, cn, thr);
  }
  if (!info == EXIT_SUCCESS) return info;

  info = cut_coordination_number(cn_max, cn, dcndr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  return EXIT_SUCCESS;
};


int get_ncoord_d4(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  TVector<double>& cn,
  TMatrix<double>& dcndr,
  bool lgrad,
  double thr
) {
  if (lgrad) {
    return dncoord_d4(mol, dist, cn, dcndr, thr);
  } 
  return ncoord_d4(mol, dist, cn, thr);
};

int ncoord_d4(const TMolecule& mol, const TMatrix<double>& dist,
              TVector<double>& cn, double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double den = 0.0;
  double countf = 0.0;
  int izp, jzp;

  for (int i = 0; i != mol.NAtoms; i++) {
    izp = mol.at(i);
    for (int j = 0; j != i; j++) {
      jzp = mol.at(j);

      r = dist(i, j);
      rcovij = rad[izp] + rad[jzp];
      rr = r / rcovij;
      den = k4 * std::exp(-pow((fabs(en[izp] - en[jzp]) + k5), 2) / k6);
      countf = den * erf_count(kn, rr);
      
      cn(i) += countf;
      cn(j) += countf;
    }
  }
  return EXIT_SUCCESS;
}


int dncoord_d4(const TMolecule& mol, const TMatrix<double>& dist,
               TVector<double>& cn, TMatrix<double>& dcndr, double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0;
  double countf = 0.0, dcountf = 0.0, den = 0.0;
  int izp, jzp;

  for (int i = 0; i != mol.NAtoms; i++) {
    izp = mol.at(i);
    for (int j = 0; j != i; j++) {
      jzp = mol.at(j);

      r = dist(i, j);
      rx = (mol.xyz(j, 0) - mol.xyz(i, 0)) / r;
      ry = (mol.xyz(j, 1) - mol.xyz(i, 1)) / r;
      rz = (mol.xyz(j, 2) - mol.xyz(i, 2)) / r;

      rcovij = rad[izp] + rad[jzp];
      rr = r / rcovij;
      den = k4 * std::exp(-pow((fabs(en[izp] - en[jzp]) + k5), 2) / k6);
      countf = den * erf_count(kn, rr);
      cn(i) += countf;
      cn(j) += countf;

      dcountf = den * derf_count(kn, rr) / rcovij;
      dcndr(j, 3 * j) += dcountf * rx;
      dcndr(j, 3 * j + 1) += dcountf * ry;
      dcndr(j, 3 * j + 2) += dcountf * rz;
      dcndr(j, 3 * i) = dcountf * rx;
      dcndr(j, 3 * i + 1) = dcountf * ry;
      dcndr(j, 3 * i + 2) = dcountf * rz;
      dcndr(i, 3 * j) = -dcountf * rx;
      dcndr(i, 3 * j + 1) = -dcountf * ry;
      dcndr(i, 3 * j + 2) = -dcountf * rz;
      dcndr(i, 3 * i) += -dcountf * rx;
      dcndr(i, 3 * i + 1) += -dcountf * ry;
      dcndr(i, 3 * i + 2) += -dcountf * rz;
    }
  }
  return EXIT_SUCCESS;
}


int ncoord_erf(const TMolecule& mol, const TMatrix<double>& dist,
               TVector<double>& cn, double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double countf = 0.0;

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != i; j++) {
      r = dist(i, j);
      rcovij = rad[mol.at(i)] + rad[mol.at(j)];
      rr = r / rcovij;
      countf = erf_count(kn, rr);
      cn(i) += countf;
      cn(j) += countf;
    }
  }

  return EXIT_SUCCESS;
}


int dncoord_erf(const TMolecule& mol, const TMatrix<double>& dist,
                TVector<double>& cn, TMatrix<double>& dcndr, double thr) {
  double r = 0.0, rcovij = 0.0, rr = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0;
  double countf = 0.0, dcountf = 0.0;

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != i; j++) {
      r = dist(i, j);
      rx = (mol.xyz(j, 0) - mol.xyz(i, 0)) / r;
      ry = (mol.xyz(j, 1) - mol.xyz(i, 1)) / r;
      rz = (mol.xyz(j, 2) - mol.xyz(i, 2)) / r;

      rcovij = rad[mol.at(i)] + rad[mol.at(j)];
      rr = r / rcovij;

      countf = erf_count(kn, rr);
      cn(i) += countf;
      cn(j) += countf;

      dcountf = derf_count(kn, rr) / rcovij;
      dcndr(j, 3 * j) += dcountf * rx;
      dcndr(j, 3 * j + 1) += dcountf * ry;
      dcndr(j, 3 * j + 2) += dcountf * rz;
      dcndr(j, 3 * i) = dcountf * rx;
      dcndr(j, 3 * i + 1) = dcountf * ry;
      dcndr(j, 3 * i + 2) = dcountf * rz;
      dcndr(i, 3 * j) = -dcountf * rx;
      dcndr(i, 3 * j + 1) = -dcountf * ry;
      dcndr(i, 3 * j + 2) = -dcountf * rz;
      dcndr(i, 3 * i) += -dcountf * rx;
      dcndr(i, 3 * i + 1) += -dcountf * ry;
      dcndr(i, 3 * i + 2) += -dcountf * rz;
    }
  }

  return EXIT_SUCCESS;
}


double erf_count(double k, double rr) {
  return 0.5 * (1.0 + std::erf(-k * (rr - 1.0)));
}


double derf_count(double k, double rr) {
  return -hlfosqrtpi * k * std::exp(-pow(k * (rr - 1.0), 2));
}


int cut_coordination_number(
  const double cn_max,
  TVector<double>& cn,
  TMatrix<double>& dcndr,
  bool lgrad
) {
  if (lgrad) {
    // TODO
  }

  for (int i = 0; i != cn.N; i++) {
    cn(i) = log_cn_cut(cn_max, cn(i));
  }

  return EXIT_SUCCESS;
};


inline double log_cn_cut(const double cn_max, const double cn) {
  return log(1.0 + exp(cn_max)) - log(1.0 + exp(cn_max - cn));
};

inline double dlog_cn_cut(const double cn_max, const double cn) {
  return exp(cn_max) / (exp(cn_max) + exp(cn));
};

}  // namespace dftd
