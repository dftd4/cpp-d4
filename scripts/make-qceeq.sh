#!/bin/bash

INCLUDE="../include"
SRC="../src"

dosed(){
	sed -n '/namespace dftd4 {/,/} \/\/ namespace dftd4/p' $1 | \
	 sed '/namespace dftd4 {/d; /} \/\/ namespace dftd4/d' | \
	 awk '/./,EOF' > $2
}

##################################################################

dosed "${INCLUDE}/dftd_eeq.h" eeq.txt

cat > qceeq.h << 'EOT'
/*
 * Electronegativity equilibration (EEQ) model for DFT-D4.
 * Contains only the essential parts for DFT-D4.
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 *
 * Extension for Fr, Ra and Actinides by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF121223)
 */

#ifndef QCEEQ_H
#define QCEEQ_H

#include "qcinpdat.h" // TGeomInput class
#include "qcmat1.h"   // TRVector and TRMatrix class

namespace dftd4 {

EOT

cat eeq.txt >> qceeq.h
rm eeq.txt

cat >> qceeq.h << 'EOT'
} // namespace dftd4

#endif // QCEEQ_H
EOT

######

dosed "${SRC}/dftd_eeq.cpp" eeq.txt

cat > qceeq.cpp << 'EOT'
/*
 * Electronegativity equilibration (EEQ) model for DFT-D4.
 * This implementation contains only the essential parts for DFT-D4.
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 *
 * Extension for Fr, Ra and Actinides by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF121223)
 */

#include <cmath>

#include "qcinpdat.h" // TGeomInput class
#include "qclineq.h"  // Linear Algebra routines
#include "qcmat1.h"   // TRVector and TRMatrix class
#include "qcmat2.h"   // BLAS routines

// CN-related logic
#include "qcncoord.h"

// always include self
#include "qceeq.h"

// wrap everything in the dftd namespace to keep it nicely confined
namespace dftd4 {

EOT

cat eeq.txt >> qceeq.cpp
rm eeq.txt

cat >> qceeq.cpp << 'EOT'
} // namespace dftd4
EOT

##########################################################################

sed -i "s/TMolecule/TGeomInput/" *.h
sed -i "s/TMolecule/TGeomInput/" *.cpp
