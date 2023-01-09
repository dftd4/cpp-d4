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
#include <cmath>
#include <iostream>

#include "dftd_cblas.h"
#include "dftd_eeq.h"
#include "dftd_dispersion.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_model.h"
#include "dftd_ncoord.h"
#include "dftd_parameters.h"

namespace dftd4 {

int get_dispersion(
  const TMolecule &mol,
  const int charge,
  const dparam &par,
  TCutoff cutoff,
  double &energy,
  double *GRAD
) {
  // setup variables
  bool lmbd = (par.s9 != 0.0);
  bool lgrad = !!GRAD;
  int info = 0;

  // distances
  TMatrix<double> dist;
  dist.NewMat(mol.NAtoms, mol.NAtoms);
  calc_distances(mol, dist);

  TVector<double> cn;        // D4 coordination number
  TVector<double> q;         // partial charges from EEQ model
  TMatrix<double> dcndr;     // derivative of D4-CN
  TMatrix<double> dqdr;      // derivative of partial charges
  TVector<double> gradient;  // derivative of dispersion energy

  cn.NewVec(mol.NAtoms);
  q.NewVec(mol.NAtoms);
  if (lgrad) {
    dcndr.NewMat(3 * mol.NAtoms, mol.NAtoms);
    dqdr.NewMat(3 * mol.NAtoms, mol.NAtoms);
    gradient.NewVec(3 * mol.NAtoms);
  }

  // calculate partial charges from EEQ model
  info = get_charges(mol, dist, charge, cutoff.cn_eeq, q, dqdr, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // get the D4 coordination number
  info = get_ncoord_d4(mol, dist, cutoff.cn, cn, dcndr, lgrad);
  if (!info == EXIT_SUCCESS) return info;


  // maximum number of reference systems
  int mref{0};
  info = get_max_ref(mol, mref); 
  if (!info == EXIT_SUCCESS) return info;

  TMatrix<double> gwvec;
  TMatrix<double> dgwdcn;
  TMatrix<double> dgwdq;
  gwvec.NewMat(mref, mol.NAtoms);
  if (lgrad) {
    dgwdcn.NewMat(mref, mol.NAtoms);
    dgwdq.NewMat(mref, mol.NAtoms);    
  }
  info = weight_references(mol, cn, q, gwvec, dgwdcn, dgwdq, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  TMatrix<double> c6;
  TMatrix<double> dc6dcn;
  TMatrix<double> dc6dq;
  c6.NewMat(mol.NAtoms, mol.NAtoms);
  if (lgrad) {
    dc6dcn.NewMat(mol.NAtoms, mol.NAtoms);
    dc6dq.NewMat(mol.NAtoms, mol.NAtoms);
  }

  info = get_atomic_c6(mol, gwvec, dgwdcn, dgwdq, c6, dc6dcn, dc6dq, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // --------------------------
  // Two-body dispersion energy
  // --------------------------

  TVector<double> dEdcn;
  TVector<double> dEdq;
  TVector<double> energies;
  energies.NewVec(mol.NAtoms);
  if (lgrad) {
    dEdcn.NewVec(mol.NAtoms);
    dEdq.NewVec(mol.NAtoms);    
  }

  info = get_dispersion2(
    mol, dist, cutoff.disp2, par, c6, dc6dcn, dc6dq, energies, dEdcn,
    dEdq, gradient, lgrad
  );
  if (!info == EXIT_SUCCESS) return info;

  if (lgrad) {
    info = BLAS_Add_Mat_x_Vec(gradient, dqdr, dEdq, false, 1.0);
    if (!info == EXIT_SUCCESS) return info;
  }

  dqdr.Delete();

  // ----------------------------
  // Three-body dispersion energy
  // ----------------------------

  if (lmbd) {
    // Three-body term is independent of charges
    for (int i = 0; i != mol.NAtoms; i++) {
      q(i) = 0.0;
    }

    // calculate weight references
    gwvec.NewMat(mref, mol.NAtoms);
    if (lgrad) {
      dgwdcn.NewMat(mref, mol.NAtoms);
      dgwdq.NewMat(mref, mol.NAtoms);    
    }
    info = weight_references(mol, cn, q, gwvec, dgwdcn, dgwdq, lgrad);
    if (!info == EXIT_SUCCESS) return info;

    cn.Delete();
    q.Delete();

    // calculate reference C6 coefficients
    c6.NewMat(mol.NAtoms, mol.NAtoms);
    if (lgrad) {
      dc6dcn.NewMat(mol.NAtoms, mol.NAtoms);
      dc6dq.NewMat(mol.NAtoms, mol.NAtoms);
    }
    info = get_atomic_c6(mol, gwvec, dgwdcn, dgwdq, c6, dc6dcn, dc6dq, lgrad);
    if (!info == EXIT_SUCCESS) return info;
    
    gwvec.Delete();
    dgwdcn.Delete();
    dgwdq.Delete();

    // calculate three-body dispersion
    info = get_dispersion3(
      mol, dist, cutoff.disp3, par, c6, dc6dcn, dc6dq, energies, dEdcn,
      dEdq, gradient, lgrad
    );
    if (!info == EXIT_SUCCESS) return info;
  } else {
    cn.Delete();
    q.Delete();
    gwvec.Delete();
    dgwdcn.Delete();
    dgwdq.Delete();
  }

  dist.Delete();
  c6.Delete();
  dc6dcn.Delete();
  dc6dq.Delete();

  if (lgrad) {
    info = BLAS_Add_Mat_x_Vec(gradient, dcndr, dEdcn, false, 1.0);
    if (!info == EXIT_SUCCESS) return info;
  }

  dcndr.Delete();
  dEdcn.Delete();
  dEdq.Delete();

  // sum up atom-wise energies
  for (int i = 0; i != mol.NAtoms; i++) {
    energy += energies(i);
  }
  energies.Delete();

  // write to input gradient
  if (lgrad) {
    for (int i = 0; i != 3*mol.NAtoms; i++) {
      GRAD[i] = gradient(i);
    }
    gradient.Delete();
  }
  
  return EXIT_SUCCESS;
}


int get_dispersion2(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient,
  bool lgrad /*= false*/
) {
  int info{0};

  if (lgrad) {
    info = get_dispersion2_derivs(
      mol, dist, cutoff, par, c6, dc6dcn, dc6dq, energy, dEdcn, dEdq, gradient
    );
  } else {
    info = get_dispersion2_energy(mol, dist, cutoff, par, c6, energy);
  }

  return info;
}

int get_dispersion2_energy(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  TVector<double>& energy
) {
  int izp{0}, jzp{0};
  double r{0.0}, r0ij{0.0}, r4r2ij{0.0}, c6ij{0.0};
  double t6{0.0}, t8{0.0}, t10{0.0};
  double e{0.0}, edisp{0.0};

  for (int iat = 0; iat != mol.NAtoms; iat++) {
    izp = mol.at(iat);
    for (int jat = 0; jat != iat; jat++) {
      jzp = mol.at(jat);
      r = dist(iat, jat);
      if (r > cutoff) continue;
      
      r4r2ij = 3.0 * r4r2[izp] * r4r2[jzp];
      r0ij = par.a1 * sqrt(r4r2ij) + par.a2;
      c6ij = c6(iat, jat);

      t6 = fdmpr_bj(6, r, r0ij);
      t8 = fdmpr_bj(8, r, r0ij);

      edisp = par.s6 * t6 + par.s8 * r4r2ij * t8;

      if (par.s10 != 0.0) {
        t10 = fdmpr_bj(10, r, r0ij);
        edisp += par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * t10;
      }

      e = -c6ij*edisp * 0.5;
      energy(iat) += e;
      energy(jat) += e;
    }
  }

  return EXIT_SUCCESS;
}

int get_dispersion2_derivs(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient
) {
  int izp{0}, jzp{0};
  double r{0.0}, r0ij{0.0}, r4r2ij{0.0}, c6ij{0.0};
  double t6{0.0}, t8{0.0}, t10{0.0};
  double dt6{0.0}, dt8{0.0}, dt10{0.0};
  double e{0.0}, edisp{0.0}, gdisp{0.0};
  double x{0.0}, y{0.0}, z{0.0};
  double dgx{0.0}, dgy{0.0}, dgz{0.0};

  for (int iat = 0; iat != mol.NAtoms; iat++) {
    izp = mol.at(iat);
    for (int jat = 0; jat != iat; jat++) {
      jzp = mol.at(jat);
      r = dist(iat, jat);
      if (r > cutoff) continue;
      
      r4r2ij = 3.0 * r4r2[izp] * r4r2[jzp];
      r0ij = par.a1 * sqrt(r4r2ij) + par.a2;
      c6ij = c6(iat, jat);

      t6 = fdmpr_bj(6, r, r0ij);
      t8 = fdmpr_bj(8, r, r0ij);

      dt6 = -6 * pow(r, 5) * pow(t6, 2);
      dt8 = -8 * pow(r, 7) * pow(t8, 2);

      edisp = par.s6 * t6 + par.s8 * r4r2ij * t8;
      gdisp = par.s6 * dt6 + par.s8 * r4r2ij * dt8;

      if (par.s10 != 0.0) {
        t10 = fdmpr_bj(10, r, r0ij);
        dt10 = -10 * pow(r, 9) * pow(t10, 2);
        edisp += par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * t10;
        gdisp += par.s10 * 49.0 / 40.0 * pow(r4r2ij, 2) * dt10;
      }

      e = -c6ij*edisp * 0.5;
      energy(iat) += e;
      energy(jat) += e;

      x = (mol.xyz(iat, 0) - mol.xyz(jat, 0)) / r;
      y = (mol.xyz(iat, 1) - mol.xyz(jat, 1)) / r;
      z = (mol.xyz(iat, 2) - mol.xyz(jat, 2)) / r;
      dgx = -c6ij * gdisp * x;
      dgy = -c6ij * gdisp * y;
      dgz = -c6ij * gdisp * z;

      dEdcn(iat) -= dc6dcn(iat, jat) * edisp;
      dEdcn(jat) -= dc6dcn(jat, iat) * edisp;

      dEdq(iat) -= dc6dq(iat, jat) * edisp;
      dEdq(jat) -= dc6dq(jat, iat) * edisp;

      gradient(3*iat) += dgx;
      gradient(3*iat + 1) += dgy;
      gradient(3*iat + 2) += dgz;
      gradient(3*jat) -= dgx;
      gradient(3*jat + 1) -= dgy;
      gradient(3*jat + 2) -= dgz;
    }
  }
  return EXIT_SUCCESS;
}


int get_dispersion3(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const dparam& par,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient,
  bool lgrad /*= false*/
) {
  return get_atm_dispersion(
    mol, dist, cutoff, par.s9, par.a1, par.a2, par.alp, c6, dc6dcn, dc6dq,
    energy, dEdcn, dEdq, gradient, lgrad
  );
}

int get_atm_dispersion(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient,
  bool lgrad /*= false*/
) {
  int info{0};

  if (lgrad) {
    info = get_atm_dispersion_derivs(
      mol, dist, cutoff, s9, a1, a2, alp, c6, dc6dcn, dc6dq, energy,
      dEdcn, dEdq, gradient
    );
  } else {
    info = get_atm_dispersion_energy(
      mol, dist, cutoff, s9, a1, a2, alp, c6, energy
    );
  }

  return info;
}

int get_atm_dispersion_energy(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double>& c6,
  TVector<double>& energy
) {
  int izp{0}, jzp{0}, kzp{0};
  double r0ij{0.0}, r0ik{0.0}, r0jk{0.0}, r0ijk{0.0};
  double c6ij{0.0}, c6ik{0.0}, c6jk{0.0}, c9ijk{0.0};
  double triple{1.0}, e{0.0};
  double rij{0.0}, rik{0.0}, rjk{0.0};
  double r2ij{0.0}, r2ik{0.0}, r2jk{0.0};
  double rijk{0.0}, r2ijk{0.0}, r3ijk{0.0};
  double ang{0.0}, fdmp{0.0};

  for (int iat = 0; iat != mol.NAtoms; iat++) {
    izp = mol.at(iat);
    for (int jat = 0; jat != iat; jat++) {
      rij = dist(iat, jat);
      if (rij > cutoff) continue; 
      r2ij = pow(rij, 2);

      jzp = mol.at(jat);
      r0ij = a1 * sqrt(3.0 * r4r2[izp] * r4r2[jzp]) + a2;
      c6ij = c6(iat, jat);

      for (int kat = 0; kat != jat; kat++) {
        rik = dist(iat, kat);
        if (rik > cutoff) continue; 
        rjk = dist(jat, kat);
        if (rjk > cutoff) continue;

        triple = triple_scale(iat, jat, kat);

        r2ik = pow(rik, 2);
        r2jk = pow(rjk, 2);

        kzp = mol.at(kat);
        r0ik = a1 * sqrt(3.0 * r4r2[izp] * r4r2[kzp]) + a2;
        r0jk = a1 * sqrt(3.0 * r4r2[jzp] * r4r2[kzp]) + a2;
        r0ijk = r0ij * r0ik * r0jk;

        c6ik = c6(iat, kat);
        c6jk = c6(jat, kat);
        c9ijk = s9 * sqrt(fabs(c6ij * c6ik * c6jk));

        rijk = rij * rik * rjk;
        r2ijk = r2ij * r2ik * r2jk;
        r3ijk = rijk * r2ijk;

        fdmp = 1.0 / (1.0 + 6.0 * pow(r0ijk / rijk, alp / 3.0));
        ang = ((0.375 * (r2ij + r2jk - r2ik) * (r2ij + r2ik - r2jk) *
               (r2ik + r2jk - r2ij) / r2ijk) + 1.0) / r3ijk;

        e = ang * fdmp * c9ijk / 3.0 * triple;
        energy(iat) += e;
        energy(jat) += e;
        energy(kat) += e;
      }
    }
  }

  return EXIT_SUCCESS;
}

double triple_scale(int ii, int jj, int kk) {
  double triple{1.0};

  if (ii == jj) {
    if (ii == kk) {
      // ii'i" -> 1/6
      triple = 1.0 / 6.0;
    } else {
      // ii'j -> 1/2
      triple = 0.5;
    }
  } else {
    if (ii != kk && jj != kk) {
      // ijk -> 1 (full)
      triple = 1.0;
    } else {
      // ijj' and iji' -> 1/2
      triple = 0.5;
    }
  }

  return triple;
}

int get_atm_dispersion_derivs(
  const TMolecule& mol,
  const TMatrix<double>& dist,
  const double cutoff,
  const double s9,
  const double a1,
  const double a2,
  const double alp,
  const TMatrix<double>& c6,
  const TMatrix<double>& dc6dcn,
  const TMatrix<double>& dc6dq,
  TVector<double>& energy,
  TVector<double>& dEdcn,
  TVector<double>& dEdq,
  TVector<double>& gradient
) {
  int izp{0}, jzp{0}, kzp{0};
  double r0ij{0.0}, r0ik{0.0}, r0jk{0.0}, r0ijk{0.0};
  double c6ij{0.0}, c6ik{0.0}, c6jk{0.0}, c9ijk{0.0};
  double triple{1.0}, e{0.0};
  double rij{0.0}, rik{0.0}, rjk{0.0};
  double xij{0.0}, xik{0.0}, xjk{0.0};
  double yij{0.0}, yik{0.0}, yjk{0.0};
  double zij{0.0}, zik{0.0}, zjk{0.0};
  double dgxij{0.0}, dgxik{0.0}, dgxjk{0.0};
  double dgyij{0.0}, dgyik{0.0}, dgyjk{0.0};
  double dgzij{0.0}, dgzik{0.0}, dgzjk{0.0};
  double r2ij{0.0}, r2ik{0.0}, r2jk{0.0};
  double rijk{0.0}, r2ijk{0.0}, r3ijk{0.0}, r5ijk{0.0};
  double ang{0.0}, dang{0.0}, fdmp{0.0}, dfdmp{0.0};
  double tmp{0.0};

  for (int iat = 0; iat != mol.NAtoms; iat++) {
    izp = mol.at(iat);
    for (int jat = 0; jat != iat; jat++) {
      rij = dist(iat, jat);
      if (rij > cutoff) continue; 
      r2ij = pow(rij, 2);

      jzp = mol.at(jat);
      r0ij = a1 * sqrt(3.0 * r4r2[izp] * r4r2[jzp]) + a2;
      c6ij = c6(iat, jat);

      for (int kat = 0; kat != jat; kat++) {
        rik = dist(iat, kat);
        if (rik > cutoff) continue; 
        rjk = dist(jat, kat);
        if (rjk > cutoff) continue;

        triple = triple_scale(iat, jat, kat);

        r2ik = pow(rik, 2);
        r2jk = pow(rjk, 2);

        kzp = mol.at(kat);
        r0ik = a1 * sqrt(3.0 * r4r2[izp] * r4r2[kzp]) + a2;
        r0jk = a1 * sqrt(3.0 * r4r2[jzp] * r4r2[kzp]) + a2;
        r0ijk = r0ij * r0ik * r0jk;

        c6ik = c6(iat, kat);
        c6jk = c6(jat, kat);
        c9ijk = -s9 * sqrt(fabs(c6ij * c6ik * c6jk));

        rijk = rij * rik * rjk;
        r2ijk = r2ij * r2ik * r2jk;
        r3ijk = rijk * r2ijk;
        r5ijk = r2ijk * r3ijk;

        tmp = pow(r0ijk / rijk, alp / 3.0);
        fdmp = 1.0 / (1.0 + 6.0 * tmp);
        ang = ((0.375 * (r2ij + r2jk - r2ik) * (r2ij + r2ik - r2jk) *
               (r2ik + r2jk - r2ij) / r2ijk) + 1.0) / r3ijk;

        e = ang * fdmp * c9ijk * triple;
        energy(iat) -= e / 3.0;
        energy(jat) -= e / 3.0;
        energy(kat) -= e / 3.0;

        // --------
        // Gradient
        // --------
        xij = (mol.xyz(jat, 0) - mol.xyz(iat, 0));
        yij = (mol.xyz(jat, 1) - mol.xyz(iat, 1));
        zij = (mol.xyz(jat, 2) - mol.xyz(iat, 2));
        xik = (mol.xyz(kat, 0) - mol.xyz(iat, 0));
        yik = (mol.xyz(kat, 1) - mol.xyz(iat, 1));
        zik = (mol.xyz(kat, 2) - mol.xyz(iat, 2));
        xjk = (mol.xyz(kat, 0) - mol.xyz(jat, 0));
        yjk = (mol.xyz(kat, 1) - mol.xyz(jat, 1));
        zjk = (mol.xyz(kat, 2) - mol.xyz(jat, 2));

        dfdmp = -2.0 * alp * tmp * pow(fdmp, 2);

        // d/drij
        dang = -0.375 * (pow(r2ij, 3) + pow(r2ij, 2) * (r2jk + r2ik)
          + r2ij * (3.0 * pow(r2jk, 2) + 2.0 * r2jk*r2ik
          + 3.0 * pow(r2ik, 2))
          - 5.0 * pow((r2jk - r2ik), 2) * (r2jk + r2ik)) / r5ijk;
        dgxij = c9ijk * (-dang*fdmp + ang*dfdmp) / r2ij * xij;
        dgyij = c9ijk * (-dang*fdmp + ang*dfdmp) / r2ij * yij;
        dgzij = c9ijk * (-dang*fdmp + ang*dfdmp) / r2ij * zij;

        // d/drik
        dang = -0.375 * (pow(r2ik, 3) + pow(r2ik, 2) * (r2jk + r2ij)
          + r2ik * (3.0 * pow(r2jk, 2) + 2.0 * r2jk * r2ij
          + 3.0 * pow(r2ij, 2))
          - 5.0 * pow((r2jk - r2ij), 2) * (r2jk + r2ij)) / r5ijk;
        dgxik = c9ijk * (-dang * fdmp + ang * dfdmp) / r2ik * xik;
        dgyik = c9ijk * (-dang * fdmp + ang * dfdmp) / r2ik * yik;
        dgzik = c9ijk * (-dang * fdmp + ang * dfdmp) / r2ik * zik;

        // d/drjk
        dang = -0.375 * (pow(r2jk, 3) + pow(r2jk, 2)*(r2ik + r2ij)
          + r2jk * (3.0 * pow(r2ik, 2) + 2.0 * r2ik * r2ij
          + 3.0 * pow(r2ij, 2))
          - 5.0 * pow((r2ik - r2ij), 2) * (r2ik + r2ij)) / r5ijk;
        dgxjk = c9ijk * (-dang * fdmp + ang * dfdmp) / r2jk * xjk;
        dgyjk = c9ijk * (-dang * fdmp + ang * dfdmp) / r2jk * yjk;
        dgzjk = c9ijk * (-dang * fdmp + ang * dfdmp) / r2jk * zjk;
        
        gradient(3*iat + 0) += - dgxij - dgxik; 
        gradient(3*iat + 1) += - dgyij - dgyik; 
        gradient(3*iat + 2) += - dgzij - dgzik;
        gradient(3*jat + 0) += dgxij - dgxjk; 
        gradient(3*jat + 1) += dgyij - dgyjk; 
        gradient(3*jat + 2) += dgzij - dgzjk; 
        gradient(3*kat + 0) += dgxik + dgxjk; 
        gradient(3*kat + 1) += dgyik + dgyjk; 
        gradient(3*kat + 2) += dgzik + dgzjk;

        dEdcn(iat) -= e * 0.5
          * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik);
        dEdcn(jat) -= e * 0.5
          * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk);
        dEdcn(kat) -= e * 0.5
          * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk);

        dEdq(iat) -= e * 0.5
          * (dc6dq(iat, jat) / c6ij + dc6dq(iat, kat) / c6ik);
        dEdq(jat) -= e * 0.5
          * (dc6dq(jat, iat) / c6ij + dc6dq(jat, kat) / c6jk);
        dEdq(kat) -= e * 0.5
          * (dc6dq(kat, iat) / c6ik + dc6dq(kat, jat) / c6jk);
      }
    }
  }

  return EXIT_SUCCESS;
}


inline double fdmpr_bj(const int n, const double r, const double c) {
  return 1.0 / (pow(r, n) + pow(c, n));
}
inline double fdmprdr_bj(const int n, const double r, const double c) {
  return -n * pow(r, n - 1) * pow(fdmpr_bj(n, r, c), 2);
}

}  // namespace dftd4
