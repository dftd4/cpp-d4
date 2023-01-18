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
#include <cmath>
#include <limits>

#include "dftd_dispersion.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_model.h"
#include "dftd_parameters.h"

namespace dftd4 {

static const double pi = 3.141592653589793238462643383279502884197;
static const double thopi = 3.0 / pi;

static const double freq[23]{0.000001, 0.050000, 0.100000, 0.200000, 0.300000,
                             0.400000, 0.500000, 0.600000, 0.700000, 0.800000,
                             0.900000, 1.000000, 1.200000, 1.400000, 1.600000,
                             1.800000, 2.000000, 2.500000, 3.000000, 4.000000,
                             5.000000, 7.500000, 10.00000};
static const double weights[23]{
  (freq[1] - freq[0]),
  (freq[1] - freq[0]) + (freq[2] - freq[1]),
  (freq[2] - freq[1]) + (freq[3] - freq[2]),
  (freq[3] - freq[2]) + (freq[4] - freq[3]),
  (freq[4] - freq[3]) + (freq[5] - freq[4]),
  (freq[5] - freq[4]) + (freq[6] - freq[5]),
  (freq[6] - freq[5]) + (freq[7] - freq[6]),
  (freq[7] - freq[6]) + (freq[8] - freq[7]),
  (freq[8] - freq[7]) + (freq[9] - freq[8]),
  (freq[9] - freq[8]) + (freq[10] - freq[9]),
  (freq[10] - freq[9]) + (freq[11] - freq[10]),
  (freq[11] - freq[10]) + (freq[12] - freq[11]),
  (freq[12] - freq[11]) + (freq[13] - freq[12]),
  (freq[13] - freq[12]) + (freq[14] - freq[13]),
  (freq[14] - freq[13]) + (freq[15] - freq[14]),
  (freq[15] - freq[14]) + (freq[16] - freq[15]),
  (freq[16] - freq[15]) + (freq[17] - freq[16]),
  (freq[17] - freq[16]) + (freq[18] - freq[17]),
  (freq[18] - freq[17]) + (freq[19] - freq[18]),
  (freq[19] - freq[18]) + (freq[20] - freq[19]),
  (freq[20] - freq[19]) + (freq[21] - freq[20]),
  (freq[21] - freq[20]) + (freq[22] - freq[21]),
  (freq[22] - freq[21])};

TD4Model::TD4Model(
  double wf_scale /*= wf_default*/,
  double ga_scale /*= ga_default*/,
  double gc_scale /*= gc_default*/
) {
  wf = wf_scale;
  ga = ga_scale;
  gc = gc_scale;
};

int TD4Model::weight_references(
  const TMolecule &mol,
  const TVector<double> &cn,
  const TVector<double> &q,
  TMatrix<double> &gwvec,
  TMatrix<double> &dgwdcn,
  TMatrix<double> &dgwdq,
  bool lgrad /*= false*/
) {
  int izp{0};
  double gw{0.0}, twf{0.0}, maxcn{0.0};
  double norm{0.0}, dnorm{0.0};
  double zi{0.0}, gi{0.0};
  double expw{0.0}, dexpw{0.0}, gwk{0.0}, dgwk{0.0};

  if (lgrad) {
    for (int iat = 0; iat != mol.NAtoms; iat++) {
      izp = mol.at(iat);
      zi = zeff[izp];
      gi = gam[izp] * gc;

      norm = 0.0;
      dnorm = 0.0;
      maxcn = 0.0;
      for (int iref = 0; iref != refn[izp]; iref++) {
        maxcn = std::max(maxcn, refcn[izp][iref]);
        for (int igw = 0; igw != refc[izp][iref]; igw++) {
          twf = (igw + 1) * wf;
          gw = weight_cn(twf, cn(iat), refcn[izp][iref]);
          norm += gw;
          dnorm += 2 * twf * (refcn[izp][iref] - cn(iat)) * gw;
        }
      }
      norm = 1.0 / norm;
      for (int iref = 0; iref != refn[izp]; iref++) {
        expw = 0.0;
        dexpw = 0.0;
        for (int igw = 0; igw != refc[izp][iref]; igw++) {
          twf = (igw + 1) * wf;
          gw = weight_cn(twf, cn(iat), refcn[izp][iref]);
          expw += gw;
          dexpw += 2 * twf * (refcn[izp][iref] - cn(iat)) * gw;
        }
        gwk = expw * norm;
        if (is_exceptional(gwk)) {
          if (refcn[izp][iref] == maxcn) {
            gwk = 1.0;
          } else {
            gwk = 0.0;
          }
        }

        gwvec(iref, iat) =
          gwk * zeta(ga, gi, refq[izp][iref] + zi, q(iat) + zi);
        dgwdq(iref, iat) =
          gwk * dzeta(ga, gi, refq[izp][iref] + zi, q(iat) + zi);

        dgwk = norm * (dexpw - expw * dnorm * norm);
        if (is_exceptional(dgwk)) { dgwk = 0.0; }
        dgwdcn(iref, iat) =
          dgwk * zeta(ga, gi, refq[izp][iref] + zi, q(iat) + zi);
      }
    }
  } else {
    for (int iat = 0; iat != mol.NAtoms; iat++) {
      izp = mol.at(iat);
      zi = zeff[izp];
      gi = gam[izp] * gc;

      norm = 0.0;
      maxcn = 0.0;
      for (int iref = 0; iref != refn[izp]; iref++) {
        maxcn = std::max(maxcn, refcn[izp][iref]);
        for (int igw = 0; igw != refc[izp][iref]; igw++) {
          twf = (igw + 1) * wf;
          norm += weight_cn(twf, cn(iat), refcn[izp][iref]);
        }
      }
      norm = 1.0 / norm;
      for (int iref = 0; iref != refn[izp]; iref++) {
        expw = 0.0;
        for (int igw = 0; igw != refc[izp][iref]; igw++) {
          twf = (igw + 1) * wf;
          expw += weight_cn(twf, cn(iat), refcn[izp][iref]);
        }
        gwk = expw * norm;
        if (std::isnan(gwk)) {
          if (refcn[izp][iref] == maxcn) {
            gwk = 1.0;
          } else {
            gwk = 0.0;
          }
        }

        gwvec(iref, iat) =
          gwk * zeta(ga, gi, refq[izp][iref] + zi, q(iat) + zi);
      }
    }
  }

  return EXIT_SUCCESS;
}

int TD4Model::get_atomic_c6(
  const TMolecule &mol,
  const TMatrix<double> &gwvec,
  const TMatrix<double> &dgwdcn,
  const TMatrix<double> &dgwdq,
  TMatrix<double> &c6,
  TMatrix<double> &dc6dcn,
  TMatrix<double> &dc6dq,
  bool lgrad /*= false*/
) {
  int izp{0}, jzp{0}, info{0};
  double refc6{0.0}, dc6{0.0};

  // maximum number of reference systems
  int mref{0};
  info = get_max_ref(mol, mref);
  if (!info == EXIT_SUCCESS) return info;

  TMatrix<double> alpha;
  alpha.NewMat(mol.NAtoms, 23 * mref);
  info = set_refalpha_eeq(mol, alpha);
  if (!info == EXIT_SUCCESS) return info;

  if (lgrad) {
    double dc6dcni{0.0}, dc6dcnj{0.0}, dc6dqi{0.0}, dc6dqj{0.0};

    for (int iat = 0; iat != mol.NAtoms; iat++) {
      izp = mol.at(iat);
      for (int jat = 0; jat != mol.NAtoms; jat++) {
        jzp = mol.at(jat);

        dc6 = 0.0;
        dc6dcni = 0.0;
        dc6dcnj = 0.0;
        dc6dqi = 0.0;
        dc6dqj = 0.0;
        for (int iref = 0; iref != refn[izp]; iref++) {
          for (int jref = 0; jref != refn[jzp]; jref++) {
            refc6 =
              thopi * trapzd(&alpha[iat][23 * iref], &alpha[jat][23 * jref]);
            dc6 += gwvec(iref, iat) * gwvec(jref, jat) * refc6;

            dc6dcni += dgwdcn(iref, iat) * gwvec(jref, jat) * refc6;
            dc6dcnj += gwvec(iref, iat) * dgwdcn(jref, jat) * refc6;
            dc6dqi += dgwdq(iref, iat) * gwvec(jref, jat) * refc6;
            dc6dqj += gwvec(iref, iat) * dgwdq(jref, jat) * refc6;
          }
        }

        c6(iat, jat) = dc6;
        c6(jat, iat) = dc6;

        dc6dcn(iat, jat) = dc6dcni;
        dc6dcn(jat, iat) = dc6dcnj;
        dc6dq(iat, jat) = dc6dqi;
        dc6dq(jat, iat) = dc6dqj;
      }
    }
  } else {
    for (int iat = 0; iat != mol.NAtoms; iat++) {
      izp = mol.at(iat);
      for (int jat = 0; jat != mol.NAtoms; jat++) {
        jzp = mol.at(jat);

        dc6 = 0.0;
        for (int iref = 0; iref != refn[izp]; iref++) {
          for (int jref = 0; jref != refn[jzp]; jref++) {
            refc6 =
              thopi * trapzd(&alpha[iat][23 * iref], &alpha[jat][23 * jref]);
            dc6 += gwvec(iref, iat) * gwvec(jref, jat) * refc6;
          }
        }

        c6(iat, jat) = dc6;
        c6(jat, iat) = dc6;
      }
    }
  }

  alpha.Delete();

  return EXIT_SUCCESS;
}

int TD4Model::set_refalpha_eeq(const TMolecule &mol, TMatrix<double> &alpha) {
  int iat{0}, is{0};
  double iz{0.0}, aiw{0.0};

  for (int i = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int ir = 0; ir != refn[iat]; ir++) {
      is = refsys[iat][ir];
      iz = zeff[is];
      for (int k = 0; k != 23; k++) {
        aiw = secscale[is] * secalpha[is][k] *
              zeta(ga, gam[is] * gc, iz, refsq[iat][ir] + iz);
        alpha(i, 23 * ir + k) = std::max(
          0.0,
          refascale[iat][ir] *
            (refalpha[iat][23 * ir + k] - refscount[iat][ir] * aiw)
        );
      }
    }
  }

  return EXIT_SUCCESS;
}

inline double trapzd(const double a[23], const double b[23]) {
  double c6 = 0.0;
  for (int w = 0; w != 23; w++)
    c6 += weights[w] * a[w] * b[w];
  return 0.5 * c6;
}

inline double weight_cn(const double wf, const double cn, const double cnref) {
  double dcn = cn - cnref;
  double arg = wf * dcn * dcn;
  return exp(-arg);
}

inline double
  zeta(const double a, const double c, const double qref, const double qmod) {
  if (qmod <= 0.0) {
    return exp(a);
  } else {
    return exp(a * (1.0 - exp(c * (1.0 - qref / qmod))));
  }
}

inline double
  dzeta(const double a, const double c, const double qref, const double qmod) {
  if (qmod <= 0.0) {
    return 0.0;
  } else {
    return -a * c * exp(c * (1.0 - qref / qmod)) * zeta(a, c, qref, qmod) *
           qref / pow(qmod, 2);
  }
}

int get_max_ref(const TMolecule &mol, int &mref) {
  mref = refn[mol.at(0)];
  for (int i = 1; i != mol.NAtoms; i++) {
    int val = refn[mol.at(i)];
    if (val > mref) mref = val;
  }

  if (mref == 0.0) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

bool is_exceptional(double val) {
  return std::isnan(val) || (fabs(val) > std::numeric_limits<double>::max());
}

} // namespace dftd4