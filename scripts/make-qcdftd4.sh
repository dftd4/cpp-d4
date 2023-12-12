#!/bin/bash

INCLUDE="../include"
SRC="../src"

cat > qcdftd4param.h << 'EOT'
/*
 * This is the DFT-D4 equivalent of copyc6, since this is an *empirical*
 * dispersion correction we can hardly avoid having this big, scary blob
 * of data lying around somewhere.
 *
 * Files like this are usually computer generated so avoid modifying it,
 * but have the appropriate tool generate it from dftd4.
 *
 * Responsible for this mess:
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 *
 * Update/Fix for parameters from periodic extension by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 *
 * Extension for Fr, Ra and Actinides by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF121223)
 */
 
#ifndef QCDFTD4PARAM_H
#define QCDFTD4PARAM_H

EOT

sed -n '/namespace/,$p' "${INCLUDE}/dftd_parameters.h" >> qcdftd4param.h

# Append the Footer:
printf "\n#endif // QCDFTD4PARAM_H\n" >> qcdftd4param.h


##################################################################
##################################################################

dosed(){
	sed -n '/namespace dftd4 {/,/} \/\/ namespace dftd4/p' $1 | \
	 sed '/namespace dftd4 {/d; /} \/\/ namespace dftd4/d' | \
	 awk '/./,EOF' > $2
}

dosed "${INCLUDE}/dftd_cutoff.h" cutoff.txt
dosed "${INCLUDE}/dftd_model.h" model.txt
dosed "${INCLUDE}/dftd_dispersion.h" dispersion.txt
dosed "${INCLUDE}/damping/dftd_rational.h" rational.txt
dosed "${INCLUDE}/damping/dftd_atm.h" atm.txt

cat > qcdftd4.h << 'EOT'
/*
 * D4(EEQ)-ATM implementation
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 */

#ifndef __QCDFTD4_H
#define __QCDFTD4_H

#include "qcinpdat.h" // TGeomInput class
#include "qcmat1.h"   // TRVector and TRMatrix class

namespace dftd4 {

EOT

######

cat >> qcdftd4.h << 'EOT'
/* --------------------------------------------------------------------------

  Cutoff (dftd_cutoff.h)
  https://github.com/dftd4/cpp-d4/blob/main/include/dftd_cutoff.h

  -------------------------------------------------------------------------- */

EOT

cat cutoff.txt >> qcdftd4.h
rm cutoff.txt

######

cat >> qcdftd4.h << 'EOT'
/* --------------------------------------------------------------------------

  D4 model (dftd_model.h)
  https://github.com/dftd4/cpp-d4/blob/main/include/dftd_model.h

  -------------------------------------------------------------------------- */

EOT

cat model.txt >> qcdftd4.h
rm model.txt

######

cat >> qcdftd4.h << 'EOT'
/* --------------------------------------------------------------------------

  Dispersion (dftd_dispersion.h)
  https://github.com/dftd4/cpp-d4/blob/main/include/dftd_dispersion.h

  -------------------------------------------------------------------------- */

EOT

cat dispersion.txt >> qcdftd4.h
rm dispersion.txt

######

cat >> qcdftd4.h << 'EOT'
/* --------------------------------------------------------------------------

  Dispersion (damping/dftd_rational.h)
  https://github.com/dftd4/cpp-d4/blob/main/include/damping/dftd_rational.h

  -------------------------------------------------------------------------- */

EOT

cat rational.txt >> qcdftd4.h
rm rational.txt

######

cat >> qcdftd4.h << 'EOT'
/* --------------------------------------------------------------------------

  Dispersion (damping/dftd_atm.h)
  https://github.com/dftd4/cpp-d4/blob/main/include/damping/dftd_atm.h

  -------------------------------------------------------------------------- */

EOT

cat atm.txt >> qcdftd4.h
rm atm.txt

######

cat >> qcdftd4.h << 'EOT'
} // namespace dftd4

/* --------------------------------------------------------------------------------------
// Calculates the EEQ charges according to the D4 paper
// https://doi.org/10.1063/1.5090222
// Bernardo de Souza, 14/09/2023
 -------------------------------------------------------------------------------------- */
int CalcEEQCharges(TRMatrix &XYZ, TIVector &ATNO, int NAtoms, int totalcharge, TRVector &q, bool printerror=true);

#endif // __QCDFTD4_H

EOT


##################################################################
##################################################################


dosed "${SRC}/dftd_cutoff.cpp" cutoff.txt
dosed "${SRC}/dftd_model.cpp" model.txt
dosed "${SRC}/dftd_dispersion.cpp" dispersion.txt
dosed "${SRC}/damping/dftd_rational.cpp" rational.txt
dosed "${SRC}/damping/dftd_atm.cpp" atm.txt

cat > qcdftd4.cpp << 'EOT'
/*
 * D4(EEQ)-ATM implementation
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 */
#include <cmath>
#include <limits>

#include "qcinpdat.h" // TGeomInput class
#include "qcmat2.h"   // BLAS routines
#include "qcmath.h"   // TVector and TMatrix class

// we cannot avoid one big, scary parameter file...
#include "qcdftd4param.h"

#include "qceeq.h"
#include "qcncoord.h"

// always include self
#include "qcdftd4.h"

/* --------------------------------------------------------------------------------------
// Calculates the EEQ charges according to the D4 paper
// https://doi.org/10.1063/1.5090222
// Bernardo de Souza, 14/09/2023
 -------------------------------------------------------------------------------------- */
int CalcEEQCharges(TRMatrix &XYZ, TIVector &ATNO, int NAtoms, int totalcharge, TRVector &q, bool printerror){

  // initialize the charges to zero
  q.Init();

  // check atomic numbers to guarantee we have all parameters
  // ghost atoms (ATNO=0) will have no charge
  for (int i=0;i<NAtoms;i++)
    if (ATNO(i)>86){
      if (printerror) printMessage("Atomic number %d detected. EEQ charges can not be calculated.\n",ATNO(i));
      return 1;
    }

  TGeomInput molecule;
  molecule.NAtoms=NAtoms;
  molecule.CC.CopyMat(XYZ);
  molecule.ATNO.CopyVec(ATNO);
  TIVector Index(NAtoms);
  for (int i=0;i<NAtoms;i++) Index(i)=i;
  TRMatrix dist;
  dist.NewMatrix(NAtoms, NAtoms);
  dftd4::calc_distances(molecule, Index, dist);

  TRMatrix dqdr;

  return dftd4::get_charges(molecule, Index, dist, totalcharge, 25, q, dqdr, false);
}

namespace dftd4 {

EOT

######

cat >> qcdftd4.cpp << 'EOT'
/* --------------------------------------------------------------------------

  Cutoff (dftd_cutoff.cpp)
  https://github.com/dftd4/cpp-d4/blob/main/src/dftd_cutoff.cpp

  -------------------------------------------------------------------------- */

EOT

cat cutoff.txt >> qcdftd4.cpp
rm cutoff.txt

######

cat >> qcdftd4.cpp << 'EOT'
/* --------------------------------------------------------------------------

  D4 model (dftd_model.cpp)
  https://github.com/dftd4/cpp-d4/blob/main/src/dftd_model.cpp

  -------------------------------------------------------------------------- */

EOT

cat model.txt >> qcdftd4.cpp
rm model.txt

######

cat >> qcdftd4.cpp << 'EOT'
/* --------------------------------------------------------------------------

  Dispersion (dftd_dispersion.cpp)
  https://github.com/dftd4/cpp-d4/blob/main/src/dftd_dispersion.cpp

  -------------------------------------------------------------------------- */

EOT

cat dispersion.txt >> qcdftd4.cpp
rm dispersion.txt

######

cat >> qcdftd4.cpp << 'EOT'
/* --------------------------------------------------------------------------

  Dispersion (damping/dftd_rational.cpp)
  https://github.com/dftd4/cpp-d4/blob/main/src/damping/dftd_rational.cpp

  -------------------------------------------------------------------------- */

EOT

cat rational.txt >> qcdftd4.cpp
rm rational.txt

######

cat >> qcdftd4.cpp << 'EOT'
/* --------------------------------------------------------------------------

  Dispersion (damping/dftd_atm.cpp)
  https://github.com/dftd4/cpp-d4/blob/main/src/damping/dftd_atm.cpp

  -------------------------------------------------------------------------- */

EOT

cat atm.txt >> qcdftd4.cpp
rm atm.txt

######

cat >> qcdftd4.cpp << 'EOT'
} // namespace dftd4
EOT


##########################################################################


sed -i "s/TMolecule/TGeomInput/" *.h
sed -i "s/TMolecule/TGeomInput/" *.cpp
