#!/bin/bash

INCLUDE="../include"
SRC="../src"

dosed(){
	sed -n '/namespace dftd4 {/,/} \/\/ namespace dftd4/p' $1 | \
	 sed '/namespace dftd4 {/d; /} \/\/ namespace dftd4/d' | \
	 awk '/./,EOF' > $2
}

##################################################################

dosed "${INCLUDE}/dftd_ncoord.h" ncoord.txt

cat > qcncoord.h << 'EOT'
/*
 * Definition of the coordination numbers used in DFT-D4 and the
 * electronegativity equilibration (EEQ) model.
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 *
 * Extension for Fr, Ra and Actinides by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF121223)
 */

#ifndef QCNCOORD_H
#define QCNCOORD_H

#include "qcinpdat.h" // TGeomInput class
#include "qcmat1.h"   // TRVector and TRMatrix class

namespace dftd4 {

EOT

cat ncoord.txt >> qcncoord.h
rm ncoord.txt

cat >> qcncoord.h << 'EOT'
} // namespace dftd4

#endif // QCNCOORD_H
EOT

######

dosed "${SRC}/dftd_ncoord.cpp" ncoord.txt

cat > qcncoord.cpp << 'EOT'
/*
 * Definition of the coordination numbers used in DFT-D4 and the
 * electronegativity equilibration (EEQ) model
 *
 * This module works on a distance matrix to avoid recalculating
 * the distances every time.
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 *
 * Extension for Fr, Ra and Actinides by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF121223)
 */

#include <cmath> // erf, exp, log, sqrt

#include "qcinpdat.h" // TGeomInput class
#include "qcmat1.h"   // TRVector and TRMatrix class

// always include self
#include "qcncoord.h"

// wrap everything in the dftd namespace to keep it nicely confined
namespace dftd4 {

EOT

cat ncoord.txt >> qcncoord.cpp
rm ncoord.txt

cat >> qcncoord.cpp << 'EOT'
} // namespace dftd4
EOT

##########################################################################

sed -i "s/TMolecule/TGeomInput/" *.h
sed -i "s/TMolecule/TGeomInput/" *.cpp
