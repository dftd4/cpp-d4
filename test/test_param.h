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
#ifndef TEST_PARAM_H
#define TEST_PARAM_H

#include <dftd_geometry.h>
#include <dftd_dispersion.h>


static const int nfuncs = 97;

static const char funcs[nfuncs][32] {
  "hf", "b-lyp", "bpbe", "b-p", "bpw", "lb94", "mpwlyp", "mpwpw", 
  "o-lyp", "opbe", "pbe", "rpbe", "revpbe", "pw86pbe", "rpw86pbe", 
  "pw91", "pwp", "x-lyp", "b97", "tpss", "revtpss", "scan", "rscan", 
  "r2scan", "b1lyp", "b3-lyp", "bh-lyp", "b1p", "b3p", "b1pw", "b3pw", 
  "o3-lyp", "revpbe0", "revpbe38", "pbe0", "pwp1", "pw1pw", "mpw1pw", 
  "mpw1lyp", "pw6b95", "tpssh", "tpss0", "x3-lyp", "m06l", "m06", 
  "wb97", "wb97x", "cam-b3lyp", "lc-blyp", "lh07tsvwn", "lh07ssvwn", 
  "lh12ctssirpw92", "lh12ctssifpw92", "lh14tcalpbe", "lh20t", 
  "b2plyp", "b2gpplyp", "mpw2plyp", "pwpb95", "dsdblyp", "dsdpbe", 
  "dsdpbeb95", "dsdpbep86", "dsdsvwn", "dodblyp", "dodpbe", "dodpbeb95", 
  "dodpbep86", "dodsvwn", "pbe02", "pbe0dh", "dftb(3ob)", "dftb(mio)", 
  "dftb(pbc)", "dftb(matsci)", "dftb(ob2)", "b1b95", "glyp", "revpbe0dh", 
  "revtpss0", "revdsd-pbep86", "revdsd-pbe", "revdsd-blyp", 
  "revdod-pbep86", "b97m", "wb97m", "pbesol", "am05", "mn12sx", 
  "hse03", "hse06", "hse12", "hse12s", "hsesol", "r2scanh", "r2scan0", 
  "r2scan50"
}; 

static const double ref_energy[nfuncs] {
  -2.8259959590382366E-01,
  -1.8404742721234266E-01,
  -1.6798937024195304E-01,
  -1.1935355217292869E-01,
  -1.7366283799329513E-01,
  -4.6430923163432714E-01,
  -1.3346580384141335E-01,
  -1.2923057292683821E-01,
  -4.2328482522988720E-01,
  -4.1963240189378620E-01,
  -8.7933519082332240E-02,
  -2.8869613111622339E-01,
  -2.5985221342803455E-01,
  -9.8893933408153903E-02,
  -9.8166674429422632E-02,
  -7.2605591514865300E-02,
  -3.3915936792186512E-02,
  -1.9922033814494208E-01,
  -1.5203207292642029E-01,
  -1.1732196176499123E-01,
  -9.4460953697244290E-02,
  -1.8508178547666513E-02,
  -3.2572027285843817E-02,
  -2.7965656650905928E-02,
  -1.4128237865557192E-01,
  -1.3614663010332351E-01,
  -9.8776916779818535E-02,
  -9.0027422431692011E-02,
  -9.4497497790429527E-02,
  -1.3584429225597902E-01,
  -1.3231891834258791E-01,
  -1.1484369432188699E-01,
  -1.8665428267713263E-01,
  -1.4700208605535828E-01,
  -7.6591931831600193E-02,
  -3.2034113023795184E-02,
  -6.6220768631983964E-02,
  -1.0533339468079540E-01,
  -1.0277805122805916E-01,
  -7.1912075732403241E-02,
  -1.0855929194374825E-01,
  -1.0118342099540853E-01,
  -1.1265405513604643E-01,
  -1.2361163873034551E-02,
  -2.0033492687049504E-02,
  -8.3715876932202084E-03,
  -1.5871986581992004E-02,
  -7.7856722692803276E-02,
  -1.2993089916236900E-02,
  -2.0792249021524978E-01,
  -2.6249584973626233E-01,
  -3.3973560197631031E-01,
  -3.7481974145720043E-01,
  -8.8553324335207614E-02,
  -5.2223863418124516E-02,
  -6.2549827081499740E-02,
  -3.9511828658811962E-02,
  -4.1587075619746859E-02,
  -6.1315118738896508E-02,
  -3.8676661542281436E-02,
  -4.7838962519922201E-02,
  -4.3721931993326975E-02,
  -1.9079922103516497E-02,
  -2.3261054412440016E-02,
  -8.5231997813884242E-02,
  -6.5803424344231073E-02,
  -5.9737404157675297E-02,
  -5.1851768992444328E-02,
  -3.1654165615043431E-02,
  -8.3257937599140481E-03,
  -4.8624385696080269E-02,
  -6.1819725698154936E-02,
  -4.5658803921709075E-02,
  -5.8629726042068023E-02,
  -7.1429204290326231E-02,
  -4.2926220854554435E-02,
  -1.0337253252219039E-01,
  -2.7118928551625038E-01,
  -1.0281292290726875E-01,
  -8.2229990650757728E-02,
  -5.6706096227474322E-02,
  -8.5751286380770605E-02,
  -8.9270144738356177E-02,
  -6.1473208957720508E-02,
  -1.2413617396073341E-01,
  -1.0557238186347086E-01,
  -3.6248687079787421E-02,
  -3.6248688696273144E-02,
  -2.4036121627224290E-02,
  -6.9745139211587009E-02,
  -7.1576872338838152E-02,
  -6.9580711257046970E-02,
  -6.2089574187874155E-02,
  -5.0422920264847002E-02,
  -2.9305044225664403E-02,
  -3.1711738228039008E-02,
  -3.4563360149465587E-02,
};


int test_dftd4_energy(
  const dftd::TMolecule &mol,
  const int &charge,
  const dftd::dparam &par,
  double &energy
);

extern int test_param(void);

#endif // TEST_PARAM_H
