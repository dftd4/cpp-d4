#!/bin/bash

INCLUDE="../include"

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

printf "\n#endif // QCDFTD4PARAM_H\n" >> qcdftd4param.h


sed -i "s/TMolecule/TGeomInput/" qcdftd4param.h
