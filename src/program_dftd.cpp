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
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "dftd_damping.h"
#include "dftd_dispersion.h"
#include "dftd_eeq.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_ncoord.h"
#include "dftd_readxyz.h"

namespace dftd {
int DFTVDW_D4(TMolecule &mol, dparam &par, int &charge, double &energy,
              double *GRAD) {
  // setup variables
  bool lverbose = false;
  bool lmbd = true;
  bool lgrad = !!GRAD;

  int info = 0;

  // this are our method constants, they could be changed, but this usually
  // breaks things
  double wf = 6.0, g_a = 3.0, g_c = 2.0;

  // local
  int ndim = 0;
  double es = 0.0;           // electrostatic energy
  TVector<double> cn;        // erf-CN for EEQ model
  TMatrix<double> dcndr;     // derivative of erf-CN
  TVector<double> q;         // partial charges from EEQ model
  TMatrix<double> dqdr;      // derivative of partial charges
  TMatrix<double> ges;       // derivative of electrostatic energy
  TMatrix<double> gradient;  // derivative of dispersion energy
  TVector<double> covcn;     // covalent CN for D4 calculation
  TMatrix<double> dcovcndr;  // derivative of covalent D4
  TVector<double> gweights;  // Gaussian weights for C6 interpolation
  TMatrix<double> c6ref;     // reference C6 coefficients
  TMatrix<double> numg;      // derivative of dispersion energy
  TMatrix<double> dist;      // distances
  dist.New(mol.NAtoms, mol.NAtoms);

  // print some information to the user
  //   if (!lsilent)
  //   {
  //      dftd4_header(lverbose);
  //      if (lverbose) gpl_license();
  //      dftd4_citation();
  //   }

  calc_distances(mol, dist);

  // calculation dimension of D4
  ndim = d4dim(mol);

  // get memory
  cn.New(mol.NAtoms);
  dcndr.New(mol.NAtoms, 3 * mol.NAtoms);
  q.New(mol.NAtoms + 1);
  dqdr.New(mol.NAtoms + 1, 3 * mol.NAtoms);
  gradient.New(mol.NAtoms, 3);
  ges.New(mol.NAtoms, 3);
  covcn.New(mol.NAtoms);
  dcovcndr.New(mol.NAtoms, 3 * mol.NAtoms);
  gweights.New(ndim);
  c6ref.New(ndim, ndim);

  // calculate partial charges by EEQ model
  info = dncoord_erf(mol, dist, cn, dcndr);
  if (!info == EXIT_SUCCESS) return info;
  info = eeq_chrgeq(mol, charge, dist, cn, dcndr, q, dqdr, es, ges, lverbose,
                    false, lgrad);
  if (!info == EXIT_SUCCESS) return info;

  // get the D4 C6 coefficients
  info = dncoord_d4(mol, dist, covcn, dcovcndr);
  if (!info == EXIT_SUCCESS) return info;
  info = d4(mol, ndim, wf, g_a, g_c, covcn, gweights, c6ref);
  if (!info == EXIT_SUCCESS) return info;

  if (!lgrad) {
    info =
        edisp(mol, dist, ndim, q, par, g_a, g_c, gweights, c6ref, lmbd, energy);
    if (!info == EXIT_SUCCESS) return info;
  } else {
    info = dispgrad(mol, dist, ndim, q, dqdr, covcn, dcovcndr, par, wf, g_a,
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
  }

  cn.Delete();
  dcndr.Delete();
  q.Delete();
  dqdr.Delete();
  gradient.Delete();
  ges.Delete();
  covcn.Delete();
  dcovcndr.Delete();
  gweights.Delete();
  c6ref.Delete();
  dist.Delete();

  return EXIT_SUCCESS;
}
};  // namespace dftd

class argparser {
 public:
  argparser(int &argc, char **argv) {
    for (int i = 1; i != argc; i++) {
      this->args.push_back(std::string(argv[i]));
    }
  }
  const std::string &getopt(const std::string &opt) const {
    std::vector<std::string>::const_iterator iter;
    iter = find(this->args.begin(), this->args.end(), opt);
    if (iter != this->args.end() && ++iter != this->args.end()) {
      return *iter;
    }
    static const std::string empty("");
    return empty;
  }
  bool getflag(const std::string &opt) const {
    return find(this->args.begin(), this->args.end(), opt) != this->args.end();
  }

 private:
  std::vector<std::string> args;
};

void dftd4_citation() {
  std::cout << R"(   Please cite this work as:
   E. Caldeweyher, C. Bannwarth and S. Grimme, J. Chem. Phys., 2017,
   147, 034112. DOI: 10.1063/1.4993215
   and
   E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher,
   C. Bannwarth and S. Grimme, J. Chem Phys, 2019, 150, 154122.
   DOI: 10.1063/1.5090222

   For a general overview on dispersion corrected mean-field methods
   we recommend:
   S. Grimme, A. Hansen, J. G. Brandenburg, C. Bannwarth, Chem. Rev.,
   2016, 116, 5105−5154.

)";
}

void gpl_license() {
  std::cout << R"(   Copyright (C) 2019-2020 S. Ehlert

   This program is free software: you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   as published by the Free Software Foundation, either version 3 of
   the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this program.
   If not, see <http://www.gnu.org/licenses/>.

)";
}

void help() {
  std::cout << R"(   Usage:
   cpp-d4 [options] <file>

   <file> is a valid Turbomole coordinate file (coordinates in Bohr) or
   in xmol format (coordinates in Ångström).

   Options:

   -f, --func <name>[/<basis>]
           calculate DFT-D4 dispersion for the given functional
   -v, --verbose
           be verbose
   -h, --help
           print this message

)";
}

int main(int argc, char **argv) {
  std::string func;
  bool lverbose{false};
  bool lmbd{true}, lgrad{false};
  dftd::dparam par;  // damping parameter for DFT-D4 calculation
  dftd::TMolecule mol;
  int info;
  double energy;
  int charge;

  // check for complete command line
  // we need at least the program name and an input file
  if (argc < 2) {
    // help();
    exit(EXIT_FAILURE);
  }
  // setup the argparser from the commandline
  argparser args(argc, argv);
  // check for help flag first
  if (args.getflag("-h") || args.getflag("--help")) {
    // help();
    exit(EXIT_SUCCESS);
  }
  if (args.getflag("-v") || args.getflag("--verbose")) {
    lverbose = true;
  }
  if (args.getflag("-g") || args.getflag("--grad")) {
    lgrad = true;
  }
  if (args.getflag("--func")) {
    func = args.getopt("--func");
    dftd::d4par(func, par, lmbd);
  }
  // last argument is assumed to filename since
  // dftd4 [options] <file>
  std::string fname{argv[argc - 1]};

  // readin the geometry file
  read_xyzfile(fname, mol);

  info = DFTVDW_D4(mol, par, charge, energy, nullptr);
  if (info != 0) return EXIT_FAILURE;

  std::cout << "Dispersion energy: " << energy << " Eh\n";

  mol.FreeMemory();

  return EXIT_SUCCESS;
}
