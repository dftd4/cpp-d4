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

/**
 * Electronegativity equilibration (EEQ) model for DFT-D4
 */

#include "dftd_eeq.h"

#include <cmath>

#include "dftd_cblas.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"

// wrap everything in the dftd namespace to keep it nicely confined
namespace dftd {

static const double xi[87]{
    0.0,        1.23695041, 1.26590957, 0.54341808, 0.99666991, 1.26691604,
    1.40028282, 1.55819364, 1.56866440, 1.57540015, 1.15056627, 0.55936220,
    0.72373742, 1.12910844, 1.12306840, 1.52672442, 1.40768172, 1.48154584,
    1.31062963, 0.40374140, 0.75442607, 0.76482096, 0.98457281, 0.96702598,
    1.05266584, 0.93274875, 1.04025281, 0.92738624, 1.07419210, 1.07900668,
    1.04712861, 1.15018618, 1.15388455, 1.36313743, 1.36485106, 1.39801837,
    1.18695346, 0.36273870, 0.58797255, 0.71961946, 0.96158233, 0.89585296,
    0.81360499, 1.00794665, 0.92613682, 1.09152285, 1.14907070, 1.13508911,
    1.08853785, 1.11005982, 1.12452195, 1.21642129, 1.36507125, 1.40340000,
    1.16653482, 0.34125098, 0.58884173, 0.68441115, 0.56999999, 0.56999999,
    0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999,
    0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999,
    0.87936784, 1.02761808, 0.93297476, 1.10172128, 0.97350071, 1.16695666,
    1.23997927, 1.18464453, 1.14191734, 1.12334192, 1.01485321, 1.12950808,
    1.30804834, 1.33689961, 1.27465977};
static const double gam[87]{
    0.0,         -0.35015861, 1.04121227,  0.09281243,  0.09412380,
    0.26629137,  0.19408787,  0.05317918,  0.03151644,  0.32275132,
    1.30996037,  0.24206510,  0.04147733,  0.11634126,  0.13155266,
    0.15350650,  0.15250997,  0.17523529,  0.28774450,  0.42937314,
    0.01896455,  0.07179178,  -0.01121381, -0.03093370, 0.02716319,
    -0.01843812, -0.15270393, -0.09192645, -0.13418723, -0.09861139,
    0.18338109,  0.08299615,  0.11370033,  0.19005278,  0.10980677,
    0.12327841,  0.25345554,  0.58615231,  0.16093861,  0.04548530,
    -0.02478645, 0.01909943,  0.01402541,  -0.03595279, 0.01137752,
    -0.03697213, 0.08009416,  0.02274892,  0.12801822,  -0.02078702,
    0.05284319,  0.07581190,  0.09663758,  0.09547417,  0.07803344,
    0.64913257,  0.15348654,  0.05054344,  0.11000000,  0.11000000,
    0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,
    0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,
    0.11000000,  0.11000000,  -0.02786741, 0.01057858,  -0.03892226,
    -0.04574364, -0.03874080, -0.03782372, -0.07046855, 0.09546597,
    0.21953269,  0.02522348,  0.15263050,  0.08042611,  0.01878626,
    0.08715453,  0.10500484};
static const double kappa[87]{
    0.0,         0.04916110,  0.10937243,  -0.12349591, -0.02665108,
    -0.02631658, 0.06005196,  0.09279548,  0.11689703,  0.15704746,
    0.07987901,  -0.10002962, -0.07712863, -0.02170561, -0.04964052,
    0.14250599,  0.07126660,  0.13682750,  0.14877121,  -0.10219289,
    -0.08979338, -0.08273597, -0.01754829, -0.02765460, -0.02558926,
    -0.08010286, -0.04163215, -0.09369631, -0.03774117, -0.05759708,
    0.02431998,  -0.01056270, -0.02692862, 0.07657769,  0.06561608,
    0.08006749,  0.14139200,  -0.05351029, -0.06701705, -0.07377246,
    -0.02927768, -0.03867291, -0.06929825, -0.04485293, -0.04800824,
    -0.01484022, 0.07917502,  0.06619243,  0.02434095,  -0.01505548,
    -0.03030768, 0.01418235,  0.08953411,  0.08967527,  0.07277771,
    -0.02129476, -0.06188828, -0.06568203, -0.11000000, -0.11000000,
    -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000,
    -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000,
    -0.11000000, -0.11000000, -0.03585873, -0.03132400, -0.05902379,
    -0.02827592, -0.07606260, -0.02123839, 0.03814822,  0.02146834,
    0.01580538,  -0.00894298, -0.05864876, -0.01817842, 0.07721851,
    0.07936083,  0.05849285};
static const double alp[87]{
    0.0,        0.55159092, 0.66205886, 0.90529132, 1.51710827, 2.86070364,
    1.88862966, 1.32250290, 1.23166285, 1.77503721, 1.11955204, 1.28263182,
    1.22344336, 1.70936266, 1.54075036, 1.38200579, 2.18849322, 1.36779065,
    1.27039703, 1.64466502, 1.58859404, 1.65357953, 1.50021521, 1.30104175,
    1.46301827, 1.32928147, 1.02766713, 1.02291377, 0.94343886, 1.14881311,
    1.47080755, 1.76901636, 1.98724061, 2.41244711, 2.26739524, 2.95378999,
    1.20807752, 1.65941046, 1.62733880, 1.61344972, 1.63220728, 1.60899928,
    1.43501286, 1.54559205, 1.32663678, 1.37644152, 1.36051851, 1.23395526,
    1.65734544, 1.53895240, 1.97542736, 1.97636542, 2.05432381, 3.80138135,
    1.43893803, 1.75505957, 1.59815118, 1.76401732, 1.63999999, 1.63999999,
    1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999,
    1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999,
    1.47055223, 1.81127084, 1.40189963, 1.54015481, 1.33721475, 1.57165422,
    1.04815857, 1.78342098, 2.79106396, 1.78160840, 2.47588882, 2.37670734,
    1.76613217, 2.66172302, 2.82773085};

static const double small = 1e-14;
static const double pi = 3.1415926535897932384626433832795029;
static const double sqrtpi = std::sqrt(pi);
static const double sqrt2pi = std::sqrt(2.0 / pi);
static const char _p_up = 'u', _p_lo = 'l';

int eeq_chrgeq(TMolecule& mol, int& charge, TMatrix<double>& dist,
               TVector<double>& cn, TMatrix<double>& dcndr, TVector<double>& q,
               TMatrix<double>& dqdr, double& energy, TMatrix<double>& gradient,
               bool lverbose, bool lgrad, bool lcpq) {
  double qtotal = 0.0;
  double es = 0.0;
  double gamij = 0.0, gamij2 = 0.0, arg = 0.0, arg2 = 0.0;
  double tmp = 0.0, dtmp = 0.0;
  double rx = 0.0, ry = 0.0, rz = 0.0, r = 0.0, r2 = 0.0;
  int info = 0;
  // lapack_int* ipiv;

  int m = mol.NAtoms + 1;
  int mm = m - 1;
  TMatrix<double> Amat;
  Amat.New(m, m);
  TVector<double> Xvec;
  Xvec.New(m);
  TVector<double> Xfac;
  Xfac.New(m);
  TVector<double> alpha;
  alpha.New(mol.NAtoms);

  for (int i = 0; i != mol.NAtoms; i++) {
    tmp = kappa[mol.at(i)] / std::sqrt(cn(i) + small);
    Xvec(i) = -xi[mol.at(i)] + tmp * cn(i);
    Xfac(i) = 0.5 * tmp;
    alpha(i) = pow(alp[mol.at(i)], 2);
  }
  Xvec(mm) = charge;

  for (int i = 0; i != mol.NAtoms; i++) {
    for (int j = 0; j != i; j++) {
      r = dist(i, j);
      gamij = 1.0 / std::sqrt(alpha(i) + alpha(j));
      arg = gamij * r;
      tmp = std::erf(arg) / r;
      Amat(i, j) = tmp;
      Amat(j, i) = tmp;
    }
    gamij = gam[mol.at(i)];
    Amat(i, i) = gamij + sqrt2pi / std::sqrt(alpha(i));
    Amat(i, mm) = 1.0;
    Amat(mm, i) = 1.0;
  }
  Amat(mm, mm) = 0.0;

  TMatrix<double> Ainv;
  Ainv.New(m, m);
  Ainv.CopyMat(Amat);

  Ainv.CopyMat(Amat);
  info = BLAS_InvertMatrix(Ainv);
  if (!info == EXIT_SUCCESS) return info;
  info = BLAS_Add_Mat_x_Vec(q, Ainv, Xvec, false, 1.0);

  qtotal = 0.0;
  for (int i = 0; i != mol.NAtoms; i++) {
    qtotal += q(i);
  }

  if (fabs(qtotal - charge) > 1.0e-8) {
    printf("Charge constraint error: %14.8f vs. %14d\n", qtotal, charge);
  }

  if (lverbose) {
    printf("    #   sym             EN              q            Aii\n");
    for (int i = 0; i != mol.NAtoms; i++) {
      printf("%5d %5d %14.8f %14.8f %14.8f\n", i, mol.at(i), -Xvec(i), q(i),
             Amat(i, i));
    }
  }

  TVector<double> scratch;
  scratch.New(m);
  for (int i = 0; i != m; i++) scratch(i) = -Xvec(i);

  info = BLAS_Add_Mat_x_Vec(scratch, Amat, q, false, 0.5);
  if (!info == EXIT_SUCCESS) return info;
  es = BLAS_Vec_x_Vec(scratch, q);
  scratch.Delete();
  energy += es;
  if (lverbose) printf("isotroptic electrostatic (IES) energy: %14.7f\n", es);

  int ThreeN = 3 * mol.NAtoms;
  TMatrix<double> dAmat;
  dAmat.New(m, ThreeN);
  TMatrix<double> dXvec;
  dXvec.New(m, ThreeN);
  TMatrix<double> Afac;
  Afac.New(mol.NAtoms, 3);
  for (int i = 0; i != mol.NAtoms; i++) {
    Afac(i, 0) = 0.0;
    Afac(i, 1) = 0.0;
    Afac(i, 2) = 0.0;
    dXvec(i, 3 * i) = dcndr(i, 3 * i) * Xfac(i);
    dXvec(i, 3 * i + 1) = dcndr(i, 3 * i + 1) * Xfac(i);
    dXvec(i, 3 * i + 2) = dcndr(i, 3 * i + 2) * Xfac(i);
    for (int j = 0; j != i; j++) {
      dXvec(i, 3 * j) = dcndr(j, 3 * i) * Xfac(i);
      dXvec(i, 3 * j + 1) = dcndr(j, 3 * i + 1) * Xfac(i);
      dXvec(i, 3 * j + 2) = dcndr(j, 3 * i + 2) * Xfac(i);
      dXvec(j, 3 * i) = dcndr(i, 3 * j) * Xfac(j);
      dXvec(j, 3 * i + 1) = dcndr(i, 3 * j + 1) * Xfac(j);
      dXvec(j, 3 * i + 2) = dcndr(i, 3 * j + 2) * Xfac(j);
      rx = mol.xyz(i, 0) - mol.xyz(j, 0);
      ry = mol.xyz(i, 1) - mol.xyz(j, 1);
      rz = mol.xyz(i, 2) - mol.xyz(j, 2);
      r2 = pow(dist(i, j), 2);
      gamij2 = 1.0 / (alpha(i) + alpha(j));  // squared
      gamij = std::sqrt(gamij2);
      arg2 = gamij2 * r2;
      dtmp = 2.0 * gamij * std::exp(-arg2) / (sqrtpi * r2) - Amat(j, i) / r2;
      Afac(i, 0) += dtmp * rx * q(j);
      Afac(i, 1) += dtmp * ry * q(j);
      Afac(i, 2) += dtmp * rz * q(j);
      Afac(j, 0) -= dtmp * rx * q(i);
      Afac(j, 1) -= dtmp * ry * q(i);
      Afac(j, 2) -= dtmp * rz * q(i);
      dAmat(j, 3 * i) = dtmp * rx * q(i);
      dAmat(j, 3 * i + 1) = dtmp * ry * q(i);
      dAmat(j, 3 * i + 2) = dtmp * rz * q(i);
      dAmat(i, 3 * j) = -dtmp * rx * q(j);
      dAmat(i, 3 * j + 1) = -dtmp * ry * q(j);
      dAmat(i, 3 * j + 2) = -dtmp * rz * q(j);
    }
  }

  if (lgrad) {
    BLAS_Add_Mat_x_Vec(gradient.p, dAmat.p, q.p, ThreeN, m, false, +1.0, 1.0);
    BLAS_Add_Mat_x_Vec(gradient.p, dXvec.p, q.p, ThreeN, m, false, -1.0, 1.0);
  }

  if (lcpq) {
    for (int i = 0; i != mol.NAtoms; i++) {
      dAmat(i, 3 * i) += Afac(i, 0);
      dAmat(i, 3 * i + 1) += Afac(i, 1);
      dAmat(i, 3 * i + 2) += Afac(i, 2);
    }

    info = BLAS_Add_Mat_x_Mat(dqdr, Ainv, dAmat, false, false, -1.0);
    if (!info == EXIT_SUCCESS) return info;
    info = BLAS_Add_Mat_x_Mat(dqdr, Ainv, dXvec, false, false, +1.0);
    if (!info == EXIT_SUCCESS) return info;
  }

  // free all memory
  Ainv.Delete();
  Amat.Delete();
  Xvec.Delete();
  Xfac.Delete();
  alpha.Delete();
  dAmat.Delete();
  dXvec.Delete();
  Afac.Delete();

  return EXIT_SUCCESS;
}

}  // namespace dftd
