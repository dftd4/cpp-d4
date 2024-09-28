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
 * Electronegativity equilibration (EEQ) model for DFT-D4.
 * This implementation contains only the essential parts for DFT-D4.
 */
#include <cmath>

#include "dftd_cblas.h"
#include "dftd_eeq.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_ncoord.h"

// wrap everything in the dftd namespace to keep it nicely confined
namespace dftd4 {

constexpr int MAXELEMENT = 104; // 103 + dummy

// Element-specific electronegativity for the EEQ charges.
static const double xi[MAXELEMENT]{
  +0.00000000, // dummy
  +1.23695041, +1.26590957, +0.54341808, +0.99666991, +1.26691604, +1.40028282,
  +1.55819364, +1.56866440, +1.57540015, +1.15056627, +0.55936220, +0.72373742,
  +1.12910844, +1.12306840, +1.52672442, +1.40768172, +1.48154584, +1.31062963,
  +0.40374140, +0.75442607, +0.76482096, +0.98457281, +0.96702598, +1.05266584,
  +0.93274875, +1.04025281, +0.92738624, +1.07419210, +1.07900668, +1.04712861,
  +1.15018618, +1.15388455, +1.36313743, +1.36485106, +1.39801837, +1.18695346,
  +0.36273870, +0.58797255, +0.71961946, +0.96158233, +0.89585296, +0.81360499,
  +1.00794665, +0.92613682, +1.09152285, +1.14907070, +1.13508911, +1.08853785,
  +1.11005982, +1.12452195, +1.21642129, +1.36507125, +1.40340000, +1.16653482,
  +0.34125098, +0.58884173, +0.68441115, +0.56999999, +0.56999999, +0.56999999,
  +0.56999999, +0.56999999, +0.56999999, +0.56999999, +0.56999999, +0.56999999,
  +0.56999999, +0.56999999, +0.56999999, +0.56999999, +0.56999999, +0.87936784,
  +1.02761808, +0.93297476, +1.10172128, +0.97350071, +1.16695666, +1.23997927,
  +1.18464453, +1.14191734, +1.12334192, +1.01485321, +1.12950808, +1.30804834,
  +1.33689961, +1.27465977, +1.06598299, +0.68184178, +1.04581665, +1.09888688,
  +1.07206461, +1.09821942, +1.10900303, +1.01039812, +1.00095966, +1.11003303,
  +1.16831853, +1.00887482, +1.05928842, +1.07672363, +1.11308426, +1.14340090,
  +1.13714110,
};

// Element-specific chemical hardnesses for the EEQ charges.
static const double gam[MAXELEMENT]{
  +0.00000000, // dummy
  -0.35015861, +1.04121227, +0.09281243, +0.09412380, +0.26629137, +0.19408787,
  +0.05317918, +0.03151644, +0.32275132, +1.30996037, +0.24206510, +0.04147733,
  +0.11634126, +0.13155266, +0.15350650, +0.15250997, +0.17523529, +0.28774450,
  +0.42937314, +0.01896455, +0.07179178, -0.01121381, -0.03093370, +0.02716319,
  -0.01843812, -0.15270393, -0.09192645, -0.13418723, -0.09861139, +0.18338109,
  +0.08299615, +0.11370033, +0.19005278, +0.10980677, +0.12327841, +0.25345554,
  +0.58615231, +0.16093861, +0.04548530, -0.02478645, +0.01909943, +0.01402541,
  -0.03595279, +0.01137752, -0.03697213, +0.08009416, +0.02274892, +0.12801822,
  -0.02078702, +0.05284319, +0.07581190, +0.09663758, +0.09547417, +0.07803344,
  +0.64913257, +0.15348654, +0.05054344, +0.11000000, +0.11000000, +0.11000000,
  +0.11000000, +0.11000000, +0.11000000, +0.11000000, +0.11000000, +0.11000000,
  +0.11000000, +0.11000000, +0.11000000, +0.11000000, +0.11000000, -0.02786741,
  +0.01057858, -0.03892226, -0.04574364, -0.03874080, -0.03782372, -0.07046855,
  +0.09546597, +0.21953269, +0.02522348, +0.15263050, +0.08042611, +0.01878626,
  +0.08715453, +0.10500484, +0.10034731, +0.15801991, -0.00071039, -0.00170887,
  -0.00133327, -0.00104386, -0.00094936, -0.00111390, -0.00125257, -0.00095936,
  -0.00102814, -0.00104450, -0.00112666, -0.00101529, -0.00059592, -0.00012585,
  -0.00140896,
};

// Element-specific CN scaling constant for the EEQ charges.
static const double kappa[MAXELEMENT]{
  +0.00000000, // dummy
  +0.04916110, +0.10937243, -0.12349591, -0.02665108, -0.02631658, +0.06005196,
  +0.09279548, +0.11689703, +0.15704746, +0.07987901, -0.10002962, -0.07712863,
  -0.02170561, -0.04964052, +0.14250599, +0.07126660, +0.13682750, +0.14877121,
  -0.10219289, -0.08979338, -0.08273597, -0.01754829, -0.02765460, -0.02558926,
  -0.08010286, -0.04163215, -0.09369631, -0.03774117, -0.05759708, +0.02431998,
  -0.01056270, -0.02692862, +0.07657769, +0.06561608, +0.08006749, +0.14139200,
  -0.05351029, -0.06701705, -0.07377246, -0.02927768, -0.03867291, -0.06929825,
  -0.04485293, -0.04800824, -0.01484022, +0.07917502, +0.06619243, +0.02434095,
  -0.01505548, -0.03030768, +0.01418235, +0.08953411, +0.08967527, +0.07277771,
  -0.02129476, -0.06188828, -0.06568203, -0.11000000, -0.11000000, -0.11000000,
  -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000,
  -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.03585873,
  -0.03132400, -0.05902379, -0.02827592, -0.07606260, -0.02123839, +0.03814822,
  +0.02146834, +0.01580538, -0.00894298, -0.05864876, -0.01817842, +0.07721851,
  +0.07936083, +0.05849285, +0.00013506, -0.00020631, +0.00473118, +0.01590519,
  +0.00369763, +0.00417543, +0.00706682, +0.00488679, +0.00505103, +0.00710682,
  +0.00463050, +0.00387799, +0.00296795, +0.00400648, +0.00548481, +0.01350400,
  +0.00675380,
};

static const double alp[MAXELEMENT]{
  +0.00000000, // dummy
  +0.55159092, +0.66205886, +0.90529132, +1.51710827, +2.86070364, +1.88862966,
  +1.32250290, +1.23166285, +1.77503721, +1.11955204, +1.28263182, +1.22344336,
  +1.70936266, +1.54075036, +1.38200579, +2.18849322, +1.36779065, +1.27039703,
  +1.64466502, +1.58859404, +1.65357953, +1.50021521, +1.30104175, +1.46301827,
  +1.32928147, +1.02766713, +1.02291377, +0.94343886, +1.14881311, +1.47080755,
  +1.76901636, +1.98724061, +2.41244711, +2.26739524, +2.95378999, +1.20807752,
  +1.65941046, +1.62733880, +1.61344972, +1.63220728, +1.60899928, +1.43501286,
  +1.54559205, +1.32663678, +1.37644152, +1.36051851, +1.23395526, +1.65734544,
  +1.53895240, +1.97542736, +1.97636542, +2.05432381, +3.80138135, +1.43893803,
  +1.75505957, +1.59815118, +1.76401732, +1.63999999, +1.63999999, +1.63999999,
  +1.63999999, +1.63999999, +1.63999999, +1.63999999, +1.63999999, +1.63999999,
  +1.63999999, +1.63999999, +1.63999999, +1.63999999, +1.63999999, +1.47055223,
  +1.81127084, +1.40189963, +1.54015481, +1.33721475, +1.57165422, +1.04815857,
  +1.78342098, +2.79106396, +1.78160840, +2.47588882, +2.37670734, +1.76613217,
  +2.66172302, +2.82773085, +1.04059593, +0.60550051, +1.22262145, +1.28736399,
  +1.44431317, +1.29032833, +1.41009404, +1.25501213, +1.15181468, +1.42010424,
  +1.43955530, +1.28565237, +1.35017463, +1.33011749, +1.30745135, +1.26526071,
  +1.34071499,
};

static const double small = 1e-14;
static const double pi = 3.1415926535897932384626433832795029;
static const double sqrtpi = std::sqrt(pi);
static const double sqrt2pi = std::sqrt(2.0 / pi);

int get_charges(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const int charge,
  const double cutoff,
  TVector<double> &q,
  TMatrix<double> &dqdr,
  bool lgrad
) {
  int info{0};
  bool lverbose{false};
  int nat = realIdx.Max() + 1;

  TVector<double> cn;    // EEQ cordination number
  TMatrix<double> dcndr; // Derivative of EEQ-CN

  cn.NewVec(nat);
  if (lgrad) dcndr.NewMat(nat, 3 * nat);

  // get the EEQ coordination number
  info = get_ncoord_erf(mol, realIdx, dist, cutoff, cn, dcndr, lgrad);
  if (info != EXIT_SUCCESS) return info;

  // corresponds to model%solve in Fortran
  info =
    eeq_chrgeq(mol, realIdx, dist, charge, cn, q, dcndr, dqdr, lgrad, lverbose);
  if (info != EXIT_SUCCESS) return info;

  dcndr.DelMat();
  cn.DelVec();

  return EXIT_SUCCESS;
};

int get_vrhs(
  const TMolecule &mol,
  const TIVector &realIdx,
  const int &charge,
  const TVector<double> &cn,
  TVector<double> &Xvec,
  TVector<double> &dXvec,
  bool lgrad
) {
  double tmp{0.0};
  int izp;
  int nat = realIdx.Max() + 1;

  if (lgrad) {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(cn(ii) + small);
      Xvec(ii) = -xi[izp] + tmp * cn(ii);
      dXvec(ii) = 0.5 * tmp;
    }
    dXvec(nat) = 0.0;
  } else {
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      izp = mol.ATNO(i);
      tmp = kappa[izp] / std::sqrt(cn(ii) + small);
      Xvec(ii) = -xi[izp] + tmp * cn(ii);
    }
  }

  // place charge at last index of xvec
  Xvec(nat) = charge;

  return EXIT_SUCCESS;
};

int get_amat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  TMatrix<double> &Amat
) {
  double gamij = 0.0;
  int mm = realIdx.Max() + 1;
  int izp, jzp;
  double alphai, alphaj;
  double tmp, r;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    izp = mol.ATNO(i);
    alphai = pow(alp[izp], 2);
    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      jzp = mol.ATNO(j);
      alphaj = pow(alp[jzp], 2);

      r = dist(ii, jj);
      gamij = 1.0 / std::sqrt(alphai + alphaj);
      tmp = std::erf(gamij * r) / r;
      Amat(ii, jj) = tmp;
      Amat(jj, ii) = tmp;
    }
    gamij = gam[izp];
    Amat(ii, ii) = gamij + sqrt2pi / alp[izp];
    Amat(ii, mm) = 1.0;
    Amat(mm, ii) = 1.0;
  }
  Amat(mm, mm) = 0.0;

  return EXIT_SUCCESS;
};

int get_damat_0d(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const TVector<double> &q,
  const TMatrix<double> &Amat,
  TMatrix<double> &dAmat,
  TMatrix<double> &atrace
) {
  double alphai, alphaj;
  double rx, ry, rz, r2;
  double arg, gam, dtmp;
  double dgx, dgy, dgz;

  for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
    ii = realIdx(i);
    if (ii < 0) continue;

    alphai = pow(alp[mol.ATNO(i)], 2);

    for (int j = 0, jj = 0; j != i; j++) {
      jj = realIdx(j);
      if (jj < 0) continue;

      alphaj = pow(alp[mol.ATNO(j)], 2);

      rx = mol.CC(i, 0) - mol.CC(j, 0);
      ry = mol.CC(i, 1) - mol.CC(j, 1);
      rz = mol.CC(i, 2) - mol.CC(j, 2);
      r2 = pow(dist(ii, jj), 2);

      gam = 1.0 / std::sqrt((alphai + alphaj));
      arg = gam * gam * r2;
      dtmp = 2.0 * gam * std::exp(-arg) / (sqrtpi * r2) - Amat(jj, ii) / r2;
      dgx = dtmp * rx;
      dgy = dtmp * ry;
      dgz = dtmp * rz;

      atrace(ii, 0) += dgx * q(jj);
      atrace(ii, 1) += dgy * q(jj);
      atrace(ii, 2) += dgz * q(jj);
      atrace(jj, 0) -= dgx * q(ii);
      atrace(jj, 1) -= dgy * q(ii);
      atrace(jj, 2) -= dgz * q(ii);

      dAmat(3 * ii, jj) = dgx * q(ii);
      dAmat(3 * ii + 1, jj) = dgy * q(ii);
      dAmat(3 * ii + 2, jj) = dgz * q(ii);
      dAmat(3 * jj, ii) = -dgx * q(jj);
      dAmat(3 * jj + 1, ii) = -dgy * q(jj);
      dAmat(3 * jj + 2, ii) = -dgz * q(jj);
    }
  }

  return EXIT_SUCCESS;
};

int eeq_chrgeq(
  const TMolecule &mol,
  const TIVector &realIdx,
  const TMatrix<double> &dist,
  const int &charge,
  const TVector<double> &cn,
  TVector<double> &qvec,
  TMatrix<double> &dcndr,
  TMatrix<double> &dqdr,
  bool lgrad /*= false*/,
  bool lverbose /*= false*/
) {
  double qtotal = 0.0;
  int info = 0;
  int n = realIdx.Max() + 1;
  int m = n + 1;

  TMatrix<double> Amat; // Coulomb matrix
  TVector<double> xvec; // x (chi) vector
  Amat.NewMat(m, m);
  xvec.NewVec(m);

  TVector<double> dxdcn; // Derivative of chi vector w.r.t. CN
  if (lgrad) dxdcn.NewVec(m);

  info = get_vrhs(mol, realIdx, charge, cn, xvec, dxdcn, lgrad);
  if (info != EXIT_SUCCESS) return info;

  info = get_amat_0d(mol, realIdx, dist, Amat);
  if (info != EXIT_SUCCESS) return info;

  TVector<double> vrhs;
  vrhs.NewVec(m);

  TMatrix<double> Ainv;
  Ainv.NewMat(m, m);
  Ainv.CopyMat(Amat);

  // solve: A Q = X (Eq.4) -> Q = Ainv X
  info = BLAS_InvertMatrix(Ainv);
  if (info != EXIT_SUCCESS) return info;

  // no return in ORCA
  BLAS_Add_Mat_x_Vec(vrhs, Ainv, xvec, false, 1.0);
  // if (info != EXIT_SUCCESS) return info;

  // remove charge constraint (make vector smaller by one) and total charge
  qtotal = 0.0;
  for (int iat = 0, ii = 0; iat != mol.NAtoms; iat++) {
    ii = realIdx(iat);
    if (ii < 0) continue;

    qvec(ii) = vrhs(ii);
    qtotal += qvec(ii);
  }

  // check total charge and additional printout
  if (fabs(qtotal - charge) > 1.0e-8) {
    printf(
      "DFT-D4: EEQ charge constraint error: %14.8f vs. %14d\n", qtotal, charge
    );
  }

  if (lverbose) {
    printf("    #   sym             EN              q            Aii\n");
    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;
      printf(
        "%5d %5d %14.8f %14.8f %14.8f\n",
        i,
        mol.ATNO(i),
        -xvec(ii),
        qvec(ii),
        Amat(ii, ii)
      );
    }
  }

  // Gradient (note that the corresponding gradient flag in Fortran is `cpq`)
  if (lgrad) {
    TMatrix<double> dAmat;
    dAmat.NewMat(3 * n, m);
    TMatrix<double> atrace;
    atrace.NewMat(m, 3);

    info = get_damat_0d(mol, realIdx, dist, vrhs, Amat, dAmat, atrace);
    if (info != EXIT_SUCCESS) return info;

    for (int i = 0, ii = 0; i != mol.NAtoms; i++) {
      ii = realIdx(i);
      if (ii < 0) continue;

      dAmat(3 * ii, ii) += atrace(ii, 0);
      dAmat(3 * ii + 1, ii) += atrace(ii, 1);
      dAmat(3 * ii + 2, ii) += atrace(ii, 2);

      for (int j = 0, jj = 0; j != mol.NAtoms; j++) {
        jj = realIdx(j);
        if (jj < 0) continue;

        dAmat(3 * jj, ii) -= dcndr(jj, 3 * ii) * dxdcn(ii);
        dAmat(3 * jj + 1, ii) -= dcndr(jj, 3 * ii + 1) * dxdcn(ii);
        dAmat(3 * jj + 2, ii) -= dcndr(jj, 3 * ii + 2) * dxdcn(ii);
      }
    }

    // we do not need these gradient-related matrices anymore
    atrace.DelMat();
    dxdcn.DelVec();

    // Ainv with last column removed
    TMatrix<double> A;
    A.NewMat(Ainv.rows, Ainv.cols - 1);
    for (int i = 0; i < Ainv.rows; i++) {
      for (int j = 0; j < Ainv.cols - 1; j++) {
        A(i, j) = Ainv(i, j);
      }
    }

    // no return in ORCA
    BLAS_Add_Mat_x_Mat(dqdr, dAmat, A, false, false, -1.0);
    // if (info != EXIT_SUCCESS) return info;

    dAmat.DelMat();
  }

  // free all memory
  Ainv.DelMat();
  Amat.DelMat();
  xvec.DelVec();

  return EXIT_SUCCESS;
}

} // namespace dftd4
