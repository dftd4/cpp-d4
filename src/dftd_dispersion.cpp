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
 * D4(EEQ)-ATM implementation
 */

#include "dftd_dispersion.h"

#include <cmath>
#include <iostream>

#include "dftd_eeq.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_ncoord.h"
#include "dftd_parameters.h"

namespace dftd {

static const double pi = 3.141592653589793238462643383279502884197;
static const double thopi = 3.0 / pi;

static const double freq[23]{0.000001, 0.050000, 0.100000, 0.200000, 0.300000,
                             0.400000, 0.500000, 0.600000, 0.700000, 0.800000,
                             0.900000, 1.000000, 1.200000, 1.400000, 1.600000,
                             1.800000, 2.000000, 2.500000, 3.000000, 4.000000,
                             5.000000, 7.500000, 10.00000};
static const double weights[23]{(freq[1] - freq[0]),
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

inline double trapzd(const double a[23], const double b[23]) {
  double c6 = 0.0;
  for (int w = 0; w != 23; w++) c6 += weights[w] * a[w] * b[w];
  return 0.5 * c6;
}

inline double cngw(const double wf, const double cn, const double cnref) {
  double dcn = cn - cnref;
  double arg = wf * dcn * dcn;
  return exp(-arg);
}

inline double zeta(const double a, const double c, const double qref,
                   const double qmod) {
  if (qmod <= 0.0) {
    return exp(a);
  } else {
    return exp(a * (1.0 - exp(c * (1.0 - qref / qmod))));
  }
}

inline double dzeta(const double a, const double c, const double qref,
                    const double qmod) {
  if (qmod <= 0.0) {
    return 0.0;
  } else {
    return -a * c * exp(c * (1.0 - qref / qmod)) * zeta(a, c, qref, qmod) *
           qref / pow(qmod, 2);
  }
}

inline double fdmpr_bj(const int n, const double r, const double c) {
  return 1.0 / (pow(r, n) + pow(c, n));
}
inline double fdmprdr_bj(const int n, const double r, const double c) {
  return -n * pow(r, n - 1) * pow(fdmpr_bj(n, r, c), 2);
}

int d4dim(const TMolecule& mol) {
  int ndim = 0;
  for (int i = 0; i != mol.NAtoms; i++) ndim += refn[mol.at(i)];
  return ndim;
}

int d4(const TMolecule& mol, int ndim, double wf, double g_a, double g_c,
       const TVector<double>& cn, TVector<double>& gw, TMatrix<double>& c6ref) {
  int iat = 0, jat = 0;
  double norm = 0.0, twf = 0.0, c6 = 0.0, maxcn = 0.0, tmp = 0.0;

  TMatrix<double> alpha;
  alpha.New(mol.NAtoms, 23*7);
  for (int i = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int ii = 0; ii != refn[iat]; ii++) {
      int is = refsys[iat][ii];
      double iz = zeff[is];
      for (int k = 0; k != 23; k++) {
        tmp = secscale[is] * secalpha[is][k]
          * zeta(g_a, gam[is] * g_c, iz, refsq[iat][ii] + iz);
        alpha(i, 23*ii + k) = std::max(0.0, refascale[iat][ii] * (refalpha[iat][23*ii + k]
              - refscount[iat][ii] * tmp));
      }
    }
  }

  for (int i = 0, k = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    norm = 0.0;
    maxcn = 0.0;
    for (int ii = 0; ii != refn[iat]; ii++) {
      maxcn = std::max(maxcn, refcn[iat][ii]);
      for (int iii = 0; iii != refc[iat][ii]; iii++) {
        twf = (iii + 1) * wf;
        norm += cngw(twf, cn(i), refcn[iat][ii]);
      }
    }
    norm = 1.0 / norm;
    for (int ii = 0; ii != refn[iat]; ii++, k++) {
      gw(k) = 0.0;
      for (int iii = 0; iii != refc[iat][ii]; iii++) {
        twf = (iii + 1) * wf;
        gw(k) += cngw(twf, cn(i), refcn[iat][ii]) * norm;
      }
      if (std::isnan(gw(k))) {
        if (refcn[iat][ii] == maxcn) {
          gw(k) = 1.0;
        } else {
          gw(k) = 0.0;
        }
      }
      for (int j = 0, l = 0; j != i; j++) {
        jat = mol.at(j);
        for (int jj = 0; jj != refn[jat]; jj++, l++) {
          c6 = thopi * trapzd(&alpha[i][23 * ii], &alpha[j][23 * jj]);
          c6ref(l, k) = c6;
          c6ref(k, l) = c6;
        }
      }
    }
  }
  alpha.Delete();

  return EXIT_SUCCESS;
}

int edisp(const TMolecule& mol, const TMatrix<double>& dist, const dparam& par, 
          int ndim, TVector<double>& q,double g_a, double g_c,
          TVector<double>& gw, TMatrix<double>& c6ref, bool lmbd,
          double& energy) {
  int iat = 0, jat = 0, ij = 0;
  int info = 0;
  double iz = 0.0, r = 0.0, r4r2ij = 0.0, cutoff = 0.0;
  double r6 = 0.0, r8 = 0.0, r10 = 0.0;
  double c6ij_ns = 0.0, c6ij = 0.0;
  double ed = 0.0, embd = 0.0;
  int packedN = mol.NAtoms * (mol.NAtoms + 1) / 2;
  TVector<double> c6ab;
  c6ab.New(packedN);
  TVector<double> zetavec;
  zetavec.New(ndim);
  TVector<double> zerovec;
  zerovec.New(ndim);

  for (int i = 0, k = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    iz = zeff[iat];
    for (int ii = 0; ii != refn[iat]; ii++, k++) {
      zetavec(k) =
          gw(k) * zeta(g_a, gam[iat] * g_c, refq[iat][ii] + iz, q(i) + iz);
      zerovec(k) = gw(k) * zeta(g_a, gam[iat] * g_c, refq[iat][ii] + iz, iz);
    }
  }

  for (int i = 0, kk = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int j = 0, ll = 0; j != i; j++) {
      ij = j + i * (i + 1) / 2;
      jat = mol.at(j);
      r = dist(i, j);
      r4r2ij = 3.0 * r4r2[iat] * r4r2[jat];
      cutoff = par.a1 * sqrt(r4r2ij) + par.a2;
      r6 = fdmpr_bj(6, r, cutoff);
      r8 = fdmpr_bj(8, r, cutoff);
      r10 = fdmpr_bj(10, r, cutoff);
      c6ij_ns = 0.0;
      c6ij = 0.0;
      for (int ii = 0, k = kk; ii != refn[iat]; ii++, k++) {
        for (int jj = 0, l = ll; jj != refn[jat]; jj++, l++) {
          c6ij_ns += zerovec(k) * zerovec(l) * c6ref(k, l);
          c6ij += zetavec(k) * zetavec(l) * c6ref(k, l);
        }
      }
      c6ab(ij) = c6ij_ns;
      ed -= c6ij * (par.s6 * r6 + par.s8 * r4r2ij * r8 +
                    par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * r10);

      ll += refn[jat];
    }
    kk += refn[iat];
  }

  if (lmbd) {
    info = apprabc(mol, dist, par, ndim, c6ab, embd);
    if (!info == 0) return info;
  }

  energy += ed + embd;

  // free all memory
  c6ab.Delete();
  zetavec.Delete();
  zerovec.Delete();

  return EXIT_SUCCESS;
}

int dispgrad(const TMolecule& mol, const TMatrix<double>& dist,
             const dparam& par, int ndim, const TVector<double>& q, 
             TMatrix<double>& dqdr, TVector<double>& cn, 
             TMatrix<double>& dcndr, double wf, double g_a, double g_c, 
             TMatrix<double>& c6ref, bool lmbd, double& energy, 
             TMatrix<double>& gradient) {
  int iat = 0, jat = 0, ij = 0;
  int info = 0;
  double norm = 0.0, dnorm = 0.0, expw = 0.0, dexpw = 0.0;
  double twf = 0.0, tgw = 0.0, gwk = 0.0, dgwk = 0.0, iz = 0.0, maxcn = 0.0;
  double disp = 0.0, c6ij = 0.0, dizij = 0.0, djzij = 0.0, dic6ij = 0.0,
         djc6ij = 0.0;
  double r = 0.0, x = 0.0, y = 0.0, z = 0.0, r4r2ij = 0.0, cutoff = 0.0,
         ed = 0.0, embd = 0.0;
  double r6 = 0.0, dr6 = 0.0, r8 = 0.0, dr8 = 0.0, r10 = 0.0, dr10 = 0.0;
  int packedN = mol.NAtoms * (mol.NAtoms + 1) / 2;
  TVector<double> zvec;
  zvec.New(ndim);
  TVector<double> gw;
  gw.New(ndim);
  TVector<double> dzvec;
  dzvec.New(ndim);
  TVector<double> dzdq;
  dzdq.New(ndim);
  TVector<double> dgw;
  dgw.New(ndim);
  TVector<double> dc6dr;
  dc6dr.New(packedN);
  TVector<double> dc6dcn;
  dc6dcn.New(mol.NAtoms);
  TVector<double> dc6dq;
  dc6dq.New(mol.NAtoms);

  for (int i = 0, k = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    iz = zeff[iat];
    norm = 0.0;
    dnorm = 0.0;
    maxcn = 0.0;
    for (int ii = 0; ii != refn[iat]; ii++) {
      maxcn = std::max(maxcn, refcn[iat][ii]);
      for (int iii = 0; iii != refc[iat][ii]; iii++) {
        twf = (iii + 1) * wf;
        tgw = cngw(twf, cn(i), refcn[iat][ii]);
        norm += tgw;
        dnorm += 2 * twf * (refcn[iat][ii] - cn(i)) * tgw;
      }
    }
    norm = 1.0 / norm;
    for (int ii = 0; ii != refn[iat]; ii++, k++) {
      dexpw = 0.0;
      expw = 0.0;
      for (int iii = 0; iii != refc[iat][ii]; iii++) {
        twf = (iii + 1) * wf;
        tgw = cngw(twf, cn(i), refcn[iat][ii]);
        expw += tgw;
        dexpw += 2 * twf * (refcn[iat][ii] - cn(i)) * tgw;
      }
      gwk = expw * norm;
      if (std::isnan(gwk)) {
        if (refcn[iat][ii] == maxcn) {
          gwk = 1.0;
        } else {
          gwk = 0.0;
        }
      }

      zvec(k) = gwk * zeta(g_a, gam[iat] * g_c, refq[iat][ii] + iz, q(i) + iz);
      gw(k) = gwk * zeta(g_a, gam[iat] * g_c, refq[iat][ii] + iz, iz);

      dgwk = dexpw * norm - expw * dnorm * pow(norm, 2);
      if (std::isnan(dgwk)) {
        dgwk = 0.0;
      }

      dzvec(k) =
          dgwk * zeta(g_a, gam[iat] * g_c, refq[iat][ii] + iz, q(i) + iz);
      dzdq(k) = gwk * dzeta(g_a, gam[iat] * g_c, refq[iat][ii] + iz, q(i) + iz);
      dgw(k) = dgwk * zeta(g_a, gam[iat] * g_c, refq[iat][ii] + iz, iz);
    }
  }

  for (int i = 0, kk = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int j = 0, ll = 0; j != i; j++) {
      ij = j + i * (i + 1) / 2;
      jat = mol.at(j);
      r = dist(i, j);

      r4r2ij = 3.0 * r4r2[iat] * r4r2[jat];
      cutoff = par.a1 * sqrt(r4r2ij) + par.a2;

      r6 = fdmpr_bj(6, r, cutoff);
      r8 = fdmpr_bj(8, r, cutoff);
      r10 = fdmpr_bj(10, r, cutoff);

      dr6 = -6 * pow(r, 5) * pow(r6, 2);
      dr8 = -8 * pow(r, 7) * pow(r8, 2);
      dr10 = -10 * pow(r, 9) * pow(r10, 2);

      c6ij = 0.0;
      dic6ij = 0.0;
      djc6ij = 0.0;
      dizij = 0.0;
      djzij = 0.0;
      for (int ii = 0, k = kk; ii != refn[iat]; ii++, k++) {
        for (int jj = 0, l = ll; jj != refn[jat]; jj++, l++) {
          c6ij += zvec(k) * zvec(l) * c6ref(k, l);
          dic6ij += dzvec(k) * zvec(l) * c6ref(k, l);
          djc6ij += zvec(k) * dzvec(l) * c6ref(k, l);
          dizij += dzdq(k) * zvec(l) * c6ref(k, l);
          djzij += zvec(k) * dzdq(l) * c6ref(k, l);
        }
      }
      disp = par.s6 * r6 + par.s8 * r4r2ij * r8 +
             par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * r10;

      // save
      dc6dq(i) += dizij * disp;
      dc6dq(j) += djzij * disp;
      dc6dcn(i) += dic6ij * disp;
      dc6dcn(j) += djc6ij * disp;
      dc6dr(ij) += c6ij * (par.s6 * dr6 + par.s8 * r4r2ij * dr8 +
                           par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * dr10);

      ed -= c6ij * disp;

      ll += refn[jat];
    }
    kk += refn[iat];
  }

  info = dabcappr(mol, dist, par, ndim, gw, dgw, c6ref, dc6dr, dc6dcn, embd);
  if (!info == 0) return info;

  energy += ed + embd;

  for (int i = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int j = 0; j != i; j++) {
      ij = j + i * (i + 1) / 2;
      r = dist(i, j);
      x = (mol.xyz(j, 0) - mol.xyz(i, 0)) / r;
      y = (mol.xyz(j, 1) - mol.xyz(i, 1)) / r;
      z = (mol.xyz(j, 2) - mol.xyz(i, 2)) / r;

      gradient(i, 0) += dc6dr(ij) * x;
      gradient(i, 1) += dc6dr(ij) * y;
      gradient(i, 2) += dc6dr(ij) * z;
      gradient(j, 0) -= dc6dr(ij) * x;
      gradient(j, 1) -= dc6dr(ij) * y;
      gradient(j, 2) -= dc6dr(ij) * z;
    }
  }

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != mol.NAtoms; j++) {
      gradient(i, 0) -= dc6dq(j) * dqdr(j, 3 * i) + dc6dcn(j) * dcndr(j, 3 * i);
      gradient(i, 1) -=
          dc6dq(j) * dqdr(j, 3 * i + 1) + dc6dcn(j) * dcndr(j, 3 * i + 1);
      gradient(i, 2) -=
          dc6dq(j) * dqdr(j, 3 * i + 2) + dc6dcn(j) * dcndr(j, 3 * i + 2);
    }
  }

  // free all memory
  zvec.Delete();
  gw.Delete();
  dzvec.Delete();
  dzdq.Delete();
  dgw.Delete();
  dc6dr.Delete();
  dc6dcn.Delete();
  dc6dq.Delete();

  return EXIT_SUCCESS;
}

int apprabc(const TMolecule& mol, const TMatrix<double>& dist,
            const dparam& par, int ndim, TVector<double>& c6ab,
            double& energy) {
  int iat = 0, jat = 0, kat = 0, ij = 0, ik = 0, jk = 0;
  double rij = 0.0, r2ij = 0.0, cij = 0.0, r4r2ij = 0.0;
  double rik = 0.0, r2ik = 0.0, cik = 0.0, r4r2ik = 0.0;
  double rjk = 0.0, r2jk = 0.0, cjk = 0.0, r4r2jk = 0.0;
  double rijk = 0.0, r2ijk = 0.0, c9ijk = 0.0, r9ijk = 0.0, r3ijk = 0.0;
  double atm = 0.0, fdmp = 0.0, cijk = 0.0, crijk = 0.0;

  for (int i = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int j = 0; j != i; j++) {
      ij = j + i * (i + 1) / 2;
      jat = mol.at(j);
      rij = dist(i, j);
      r2ij = pow(rij, 2);
      r4r2ij = 3.0 * r4r2[iat] * r4r2[jat];
      cij = par.a1 * sqrt(r4r2ij) + par.a2;
      for (int k = 0; k != j; k++) {
        ik = k + i * (i + 1) / 2;
        jk = k + j * (j + 1) / 2;
        kat = mol.at(k);
        rik = dist(i, k);
        r2ik = pow(rik, 2);
        r4r2ik = 3.0 * r4r2[iat] * r4r2[kat];
        cik = par.a1 * sqrt(r4r2ik) + par.a2;
        rjk = dist(j, k);
        r2jk = pow(rjk, 2);
        r4r2jk = 3.0 * r4r2[jat] * r4r2[kat];
        cjk = par.a1 * sqrt(r4r2jk) + par.a2;
        r2ijk = r2ij * r2ik * r2jk;
        rijk = rij * rik * rjk;
        r3ijk = rijk * r2ijk;
        cijk = cij * cik * cjk;
        c9ijk = par.s9 * sqrt(std::max(0.0, c6ab(ij) * c6ab(jk) * c6ab(ik)));
        atm = (0.375 * (r2ij + r2jk - r2ik) * (r2ij + r2ik - r2jk) *
               (r2ik + r2jk - r2ij) / r2ijk) +
              1.0;
        crijk = pow(cijk / rijk, 1.0 / 3.0);
        fdmp = 1.0 / (1.0 + 6.0 * pow(crijk, par.alp));
        r9ijk = atm / r3ijk * fdmp;
        energy += c9ijk * r9ijk;
      }
    }
  }

  return EXIT_SUCCESS;
}

int dabcappr(const TMolecule& mol, const TMatrix<double>& dist, 
             const dparam& par, int ndim, TVector<double>& gw, 
             TVector<double>& dgw, TMatrix<double>& c6ref, 
             TVector<double>& dc6dr, TVector<double>& dc6dcn, double& energy) {
  int iat = 0, jat = 0, kat = 0, ij = 0, ik = 0, jk = 0;
  double c6ij = 0.0, dic6ij = 0.0, djc6ij = 0.0, dtmp = 0.0, rij = 0.0,
         r2ij = 0.0, cij = 0.0, dijfdmp = 0.0, dijatm = 0.0, r4r2ij = 0.0,
         rik = 0.0, r2ik = 0.0, cik = 0.0, dikfdmp = 0.0, dikatm = 0.0,
         r4r2ik = 0.0, rjk = 0.0, r2jk = 0.0, cjk = 0.0, djkfdmp = 0.0,
         djkatm = 0.0, r4r2jk = 0.0, rijk = 0.0, r2ijk = 0.0, c9ijk = 0.0,
         r9ijk = 0.0, atm = 0.0, fdmp = 0.0, cijk = 0.0, dic9ijk = 0.0,
         djc9ijk = 0.0, dkc9ijk = 0.0, r3ijk = 0.0, crijk = 0.0;

  TVector<double> c6ab;
  c6ab.New(mol.NAtoms * (mol.NAtoms + 1) / 2);
  TMatrix<double> dc6ab;
  dc6ab.New(mol.NAtoms, mol.NAtoms);

  for (int i = 0, kk = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int j = 0, ll = 0; j != i; j++) {
      ij = j + i * (i + 1) / 2;
      jat = mol.at(j);

      c6ij = 0.0;
      dic6ij = 0.0;
      djc6ij = 0.0;
      for (int ii = 0, k = kk; ii != refn[iat]; ii++, k++) {
        for (int jj = 0, l = ll; jj != refn[jat]; jj++, l++) {
          c6ij += gw(k) * gw(l) * c6ref(k, l);
          dic6ij += dgw(k) * gw(l) * c6ref(k, l);
          djc6ij += gw(k) * dgw(l) * c6ref(k, l);
        }
      }
      // save
      c6ab(ij) = c6ij;
      dc6ab(i, j) = dic6ij;
      dc6ab(j, i) = djc6ij;

      ll += refn[jat];
    }
    kk += refn[iat];
  }

  for (int i = 0; i != mol.NAtoms; i++) {
    iat = mol.at(i);
    for (int j = 0; j != i; j++) {
      ij = j + i * (i + 1) / 2;
      jat = mol.at(j);
      rij = dist(i, j);
      r2ij = pow(rij, 2);
      r4r2ij = 3.0 * r4r2[iat] * r4r2[jat];
      cij = par.a1 * sqrt(r4r2ij) + par.a2;
      for (int k = 0; k != j; k++) {
        ik = k + i * (i + 1) / 2;
        jk = k + j * (j + 1) / 2;
        kat = mol.at(k);
        rik = dist(i, k);
        r2ik = pow(rik, 2);
        r4r2ik = 3.0 * r4r2[iat] * r4r2[kat];
        cik = par.a1 * sqrt(r4r2ik) + par.a2;
        rjk = dist(j, k);
        r2jk = pow(rjk, 2);
        r4r2jk = 3.0 * r4r2[jat] * r4r2[kat];
        cjk = par.a1 * sqrt(r4r2jk) + par.a2;
        r2ijk = r2ij * r2ik * r2jk;
        rijk = rij * rik * rjk;
        r3ijk = rijk * r2ijk;
        cijk = cij * cik * cjk;

        atm = ((0.375 * (r2ij + r2jk - r2ik) * (r2ij + r2ik - r2jk) *
                (r2ik + r2jk - r2ij) / r2ijk) +
               1.0) /
              r3ijk;
        dijatm = -0.375 *
                 (pow(r2ij, 3) + pow(r2ij, 2) * (r2jk + r2ik) +
                  r2ij * (3.0 * pow(r2jk, 2) + 2.0 * r2jk * r2ik +
                          3.0 * pow(r2ik, 2)) -
                  5.0 * pow(r2jk - r2ik, 2) * (r2jk + r2ik)) /
                 (r2ijk * r3ijk) / rij;
        djkatm = -0.375 *
                 (pow(r2jk, 3) + pow(r2jk, 2) * (r2ik + r2ij) +
                  r2jk * (3.0 * pow(r2ik, 2) + 2.0 * r2ik * r2ij +
                          3.0 * pow(r2ij, 2)) -
                  5.0 * pow(r2ik - r2ij, 2) * (r2ik + r2ij)) /
                 (r2ijk * r3ijk) / rjk;
        dikatm = -0.375 *
                 (pow(r2ik, 3) + pow(r2ik, 2) * (r2jk + r2ij) +
                  r2ik * (3.0 * pow(r2jk, 2) + 2.0 * r2jk * r2ij +
                          3.0 * pow(r2ij, 2)) -
                  5.0 * pow(r2jk - r2ij, 2) * (r2jk + r2ij)) /
                 (r2ijk * r3ijk) / rik;

        crijk = pow(cijk / rijk, 1.0 / 3.0);
        fdmp = 1.0 / (1.0 + 6.0 * pow(crijk, par.alp));
        dtmp = -pow(fdmp, 2) * 6.0 / 3.0 * par.alp * pow(crijk, par.alp);
        dijfdmp = dtmp / rij;
        djkfdmp = dtmp / rjk;
        dikfdmp = dtmp / rik;

        c9ijk = par.s9 * sqrt(std::max(0.0, c6ab(ij) * c6ab(jk) * c6ab(ik)));
        r9ijk = atm * fdmp;

        energy += c9ijk * r9ijk;

        dc6dr(ij) += (atm * dijfdmp - dijatm * fdmp) * c9ijk;
        dc6dr(ik) += (atm * dikfdmp - dikatm * fdmp) * c9ijk;
        dc6dr(jk) += (atm * djkfdmp - djkatm * fdmp) * c9ijk;
        dic9ijk = dc6ab(i, j) / c6ab(ij) + dc6ab(i, k) / c6ab(ik);
        djc9ijk = dc6ab(j, i) / c6ab(ij) + dc6ab(j, k) / c6ab(jk);
        dkc9ijk = dc6ab(k, j) / c6ab(jk) + dc6ab(k, i) / c6ab(ik);
        dc6dcn(i) -= 0.5 * c9ijk * r9ijk * dic9ijk;
        dc6dcn(j) -= 0.5 * c9ijk * r9ijk * djc9ijk;
        dc6dcn(k) -= 0.5 * c9ijk * r9ijk * dkc9ijk;
      }
    }
  }

  // free all memory
  c6ab.Delete();
  dc6ab.Delete();

  return EXIT_SUCCESS;
}

int DFTVDW_D4(const TMolecule &mol, const dparam &par, const int &charge,
              double &energy, double *GRAD) {
  // setup variables
  bool lverbose = false;
  bool lmbd = true;
  bool lgrad = !!GRAD;

  int info = 0;
  int ndim = 0;

  // this are our method constants, they could be changed, but this usually
  // breaks things
  double wf = 6.0, g_a = 3.0, g_c = 2.0;

  // distances
  TMatrix<double> dist;
  dist.New(mol.NAtoms, mol.NAtoms);
  calc_distances(mol, dist);

  // calculation dimension of D4
  ndim = d4dim(mol);
  

  double es = 0.0;           // electrostatic energy
  TVector<double> covcn;     // D4 coordination number
  TVector<double> cn;        // EEQ cordination number
  TVector<double> q;         // partial charges from EEQ model
  TVector<double> gweights;  // Gaussian weights for C6 interpolation
  TMatrix<double> c6ref;     // reference C6 coefficients
  TMatrix<double> numg;      // derivative of dispersion energy

  // get memory
  c6ref.New(ndim, ndim);
  covcn.New(mol.NAtoms);
  cn.New(mol.NAtoms);
  gweights.New(ndim);
  q.New(mol.NAtoms + 1);


  TMatrix<double> dcndr;     // derivative of erf-CN
  TMatrix<double> dcovcndr;  // derivative of covalent D4
  TMatrix<double> dqdr;      // derivative of partial charges
  TMatrix<double> ges;       // derivative of electrostatic energy
  TMatrix<double> gradient;  // derivative of dispersion energy
  if (lgrad) {
    dcndr.New(mol.NAtoms, 3 * mol.NAtoms);
    dcovcndr.New(mol.NAtoms, 3 * mol.NAtoms);
    dqdr.New(mol.NAtoms + 1, 3 * mol.NAtoms);
    ges.New(mol.NAtoms, 3);
    gradient.New(mol.NAtoms, 3);
  } 
  
  // get the EEQ coordination number
  info = get_ncoord_erf(mol, dist, cn, dcndr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // calculate partial charges by EEQ model
  info = eeq_chrgeq(
    mol, dist, charge, cn, q, es, dcndr, dqdr, ges, lgrad, lverbose, false
  );
  if (!info == EXIT_SUCCESS) return info;

  // get the D4 coordination number
  info = get_ncoord_d4(mol, dist, covcn, dcovcndr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // D4 weights and c6 references
  info = d4(mol, ndim, wf, g_a, g_c, covcn, gweights, c6ref);
  if (!info == EXIT_SUCCESS) return info;

  if (!lgrad) {
    info =
        edisp(mol, dist, par, ndim, q, g_a, g_c, gweights, c6ref, lmbd, energy);
    if (!info == EXIT_SUCCESS) return info;
  } else {
    info = dispgrad(mol, dist, par, ndim, q, dqdr, covcn, dcovcndr, wf, g_a,
                    g_c, c6ref, lmbd, energy, gradient);
    if (!info == EXIT_SUCCESS) return info;
    // add to gradient
    for (int i = 0, ij = 0; i != mol.NAtoms; i++) {
      for (int j = 0; j != 3; j++, ij++) {
        // printf("%14.8f", gradient(i,j));
        GRAD[ij] += gradient(i, j);
      }
      // printf("\n");
    }

    dcndr.Delete();
    dcovcndr.Delete();
    dqdr.Delete();
    ges.Delete();
    gradient.Delete();
  }

  c6ref.Delete();
  cn.Delete();
  covcn.Delete();
  dist.Delete();
  gweights.Delete();
  q.Delete();

  return EXIT_SUCCESS;
}

}  // namespace dftd
