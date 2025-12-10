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
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "dftd_cutoff.h"
#include "dftd_damping.h"
#include "dftd_dispersion.h"
#include "dftd_geometry.h"
#include "dftd_matrix.h"
#include "dftd_model.h"
#include "dftd_readxyz.h"

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
      return find(this->args.begin(), this->args.end(), opt) !=
             this->args.end();
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
  std::cout << R"(   Copyright (C) 2019 S. Ehlert, M. Friede

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

   <file> is a valid xmol file (coordinates in Ångström).
   (Turbomole coordinate file (coordinates in Bohr) not yet supported!)


   Options:

   -f, --func <name>[/<basis>]
           calculate DFT-D4 dispersion for the given functional
   -v, --verbose
           be verbose
   -h, --help
           print this message
   -g, --grad
           calculate dispersion gradient
)";
}

int main(int argc, char **argv) {
  std::string func;
  bool lmbd{true}, lgrad{false};
  dftd4::dparam par; // damping parameter for DFT-D4 calculation
  dftd4::TMolecule mol;
  int info{0};
  int charge{0};
  double energy{0.0};

  // check for complete command line
  // we need at least the program name and an input file
  if (argc < 2) {
    help();
    exit(EXIT_FAILURE);
  }

  // setup the argparser from the commandline
  argparser args(argc, argv);

  // check for help flag first
  if (args.getflag("-h") || args.getflag("--help")) {
    help();
    exit(EXIT_SUCCESS);
  }

  // get other flags
  // bool lverbose{false};
  // if (args.getflag("-v") || args.getflag("--verbose")) { lverbose = true; }
  if (args.getflag("-g") || args.getflag("--grad")) { lgrad = true; }

  if (args.getflag("--func")) {
    func = args.getopt("--func");
  } else if (args.getflag("-f")) {
    func = args.getopt("-f");
  } else {
    printf(
      "WARNING: No functional given (via '--func'). Defaulting to PBE.\n\n"
    );
    func = "pbe";
  }
  dftd4::d4par(func, par, lmbd);

  // last argument is assumed to filename since
  // dftd4 [options] <file>
  std::string fname{argv[argc - 1]};

  // readin the geometry file
  read_xyzfile(fname, mol);

  // initialize default cutoffs and default D4 model
  dftd4::TCutoff cutoff;
  dftd4::TD4Model d4;

  // analytical gradient
  double *d4grad;
  if (lgrad) {
    d4grad = new double[3 * mol.NAtoms];
    for (int i = 0; i < 3 * mol.NAtoms; i++) {
      d4grad[i] = 0.0;
    }
  } else {
    d4grad = nullptr;
  }

  info = dftd4::get_dispersion(mol, charge, d4, par, cutoff, energy, d4grad);
  if (info != EXIT_SUCCESS) return info;

  // Print results
  printf("Dispersion energy: %.15e Eh\n", energy);

  if (lgrad) {
    printf("\nDispersion gradient\n");
    int count{0};
    for (int i = 0; i < mol.NAtoms; i++) {
      printf(
        "%3i  :  %+14.9le  %+14.9le  %+14.9le\n",
        i + 1,
        d4grad[count],
        d4grad[count + 1],
        d4grad[count + 2]
      );
      count = count + 3;
    };
  }

  mol.FreeMemory();
  delete[] d4grad;

  return EXIT_SUCCESS;
}
