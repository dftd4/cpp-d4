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
#include <map>
#include <string>

#include "dftd_damping.h"
#include "dftd_dispersion.h"

namespace dftd {

enum dfunc {
  none = 0,
  hf = 1,
  blyp = 2,
  bpbe = 3,
  bp = 4,
  bpw = 5,
  lb94 = 6,
  mpwlyp = 7,
  mpwpw = 8,
  olyp = 9,
  opbe = 10,
  pbe = 11,
  rpbe = 12,
  revpbe = 13,
  pw86pbe = 14,
  rpw86pbe = 15,
  pw91 = 16,
  pwp = 17,
  xlyp = 18,
  b97 = 19,
  tpss = 20,
  revtpss = 21,
  scan = 22,
  b1lyp = 23,
  b3lyp = 24,
  bhlyp = 25,
  b1p = 26,
  b3p = 27,
  b1pw = 28,
  b3pw = 29,
  o3lyp = 30,
  revpbe0 = 31,
  revpbe38 = 32,
  pbe0 = 33,
  pwp1 = 34,
  pw1pw = 35,
  mpw1pw = 36,
  mpw1lyp = 37,
  pw6b95 = 38,
  tpssh = 39,
  tpss0 = 40,
  x3lyp = 41,
  m06l = 42,
  m06 = 43,
  m062x = 44,
  wb97 = 45,
  wb97x = 46,
  camb3lyp = 47,
  lcblyp = 48,
  lh07tsvwn = 49,
  lh07ssvwn = 50,
  lh12ctssirpw92 = 51,
  lh12ctssifpw92 = 52,
  lh14tcalpbe = 53,
  b2plyp = 54,
  b2gpplyp = 55,
  mpw2plyp = 56,
  pwpb95 = 57,
  dsdblyp = 58,
  dsdpbe = 59,
  dsdpbeb95 = 60,
  dsdpbep86 = 61,
  dsdsvwn = 62,
  dodblyp = 63,
  dodpbe = 64,
  dodpbeb95 = 65,
  dodpbep86 = 66,
  dodsvwn = 67,
  pbe0_2 = 68,
  pbe0_dh = 69,
  hf3c = 70,
  hf3cv = 71,
  pbeh3c = 72,
  b973c = 73,
  hsesol = 74,
  pwgga = 75,
  dftb3 = 76,
  hcth120 = 77,
  ptpss = 78,
  lcwpbe = 79,
  bmk = 80,
  b1b95 = 81,
  pwb6k = 82,
  otpss = 83,
  ssb = 84,
  revssb = 85,
  pbesol = 86,
  hse06 = 87,
  pbexalpha = 88,
  pbehpbe = 89,
  hcth407 = 90,
  n12 = 91,
  pkzb = 92,
  thcth = 93,
  m11l = 94,
  mn15l = 95,
  mpwb1k = 96,
  mpw1kcis = 97,
  mpwkcis1k = 98,
  pbeh1pbe = 99,
  pbe1kcis = 100,
  b97_1 = 101,
  b97_2 = 102,
  b98 = 103,
  hiss = 104,
  hse03 = 105,
  revtpssh = 106,
  tpss1kcis = 107,
  m05 = 108,
  m052x = 109,
  m08hx = 110,
  lcwhpbe = 111,
  mn12l = 112,
  tauhcthhyb = 113,
  sogga11x = 114,
  n12sx = 115,
  mn12sx = 116,
  mn15 = 117,
  glyp = 118,
  bop = 119,
  mpw1b95 = 120,
  revpbe0dh = 121,
  revtpss0 = 122,
  revdsdpbep86 = 123,
  revdsdpbe = 124,
  revdsdblyp = 125,
  revdodpbep86 = 126,
  am05 = 127,
  hse12 = 128,
  hse12s = 129,
  rscan = 130,
  r2scan = 131,
  r2scanh = 132,
  r2scan0 = 133,
  r2scan50 = 134,
  lh20t = 135,
  dftb_mio = 137,
  dftb_ob2 = 138,
  dftb_matsci = 139,
  dftb_pbc = 140,
  b97d = 141,
  b97m = 142,
  wb97m = 143,
};

dfunc get_dfunc(std::string name) {
  static const std::map<std::string, dfunc> dfuncs{
      {"hf", hf},
      {"b-lyp", blyp},
      {"blyp", blyp},
      {"bpbe", bpbe},
      {"b-p", bp},
      {"bp86", bp},
      {"bp", bp},
      {"b-p86", bp},
      {"bpw", bpw},
      {"b-pw", bpw},
      {"lb94", lb94},
      {"mpwlyp", mpwlyp},
      {"mpw-lyp", mpwlyp},
      {"mpwpw", mpwpw},
      {"mpw-pw", mpwpw},
      {"mpwpw91", mpwpw},
      {"o-lyp", olyp},
      {"olyp", olyp},
      {"opbe", opbe},
      {"pbe", pbe},
      {"rpbe", rpbe},
      {"revpbe", revpbe},
      {"pw86pbe", pw86pbe},
      {"rpw86pbe", rpw86pbe},
      {"pw91", pw91},
      {"pwp", pwp},
      {"pw-p", pwp},
      {"pw91p86", pwp},
      {"x-lyp", xlyp},
      {"xlyp", xlyp},
      {"b97", b97},
      {"tpss", tpss},
      {"revtpss", revtpss},
      {"scan", scan},
      {"b1lyp", b1lyp},
      {"b1-lyp", b1lyp},
      {"b3-lyp", b3lyp},
      {"b3lyp", b3lyp},
      {"bh-lyp", bhlyp},
      {"bhlyp", bhlyp},
      {"b1p", b1p},
      {"b1-p", b1p},
      {"b1p86", b1p},
      {"b3p", b3p},
      {"b3-p", b3p},
      {"b3p86", b3p},
      {"b1pw", b1pw},
      {"b1-pw", b1pw},
      {"b1pw91", b1pw},
      {"b3pw", b3pw},
      {"b3-pw", b3pw},
      {"b3pw91", b3pw},
      {"o3-lyp", o3lyp},
      {"o3lyp", o3lyp},
      {"revpbe0", revpbe0},
      {"revpbe38", revpbe38},
      {"pbe0", pbe0},
      {"pwp1", pwp1},
      {"pw1pw", pw1pw},
      {"pw1-pw", pw1pw},
      {"mpw1pw", mpw1pw},
      {"mpw1-pw", mpw1pw},
      {"mpw1pw91", mpw1pw},
      {"mpw1lyp", mpw1lyp},
      {"mpw1-lyp", mpw1lyp},
      {"pw6b95", pw6b95},
      {"tpssh", tpssh},
      {"tpss0", tpss0},
      {"x3-lyp", x3lyp},
      {"x3lyp", x3lyp},
      {"m06l", m06l},
      {"m06", m06},
      {"m06-2x", m062x},
      {"m062x", m062x},
      {"b97d", b97d},
      {"b97m", b97m},
      {"wb97", wb97},
      {"ωb97", wb97},
      {"omegab97", wb97},
      {"wb97x", wb97x},
      {"ωb97x", wb97x},
      {"omegab97x", wb97x},
      {"wb97m", wb97m},
      {"ωb97m", wb97m},
      {"omegab97m", wb97m},
      {"cam-b3lyp", camb3lyp},
      {"lc-blyp", lcblyp},
      {"lh07tsvwn", lh07tsvwn},
      {"lh07t-svwn", lh07tsvwn},
      {"lh07ssvwn", lh07ssvwn},
      {"lh07s-svwn", lh07ssvwn},
      {"lh12ctssirpw92", lh12ctssirpw92},
      {"lh12ct-ssirpw92", lh12ctssirpw92},
      {"lh12ctssifpw92", lh12ctssifpw92},
      {"lh12ct-ssifpw92", lh12ctssifpw92},
      {"lh14tcalpbe", lh14tcalpbe},
      {"lh14t-calpbe", lh14tcalpbe},
      {"lh20t", lh20t},
      {"b2plyp", b2plyp},
      {"b2-plyp", b2plyp},
      {"b2gpplyp", b2gpplyp},
      {"b2gp-plyp", b2gpplyp},
      {"mpw2plyp", mpw2plyp},
      {"pwpb95", pwpb95},
      {"dsdblyp", dsdblyp},
      {"dsd-blyp", dsdblyp},
      {"dsdpbe", dsdpbe},
      {"dsd-pbe", dsdpbe},
      {"dsdpbeb95", dsdpbeb95},
      {"dsd-pbeb95", dsdpbeb95},
      {"dsdpbep86", dsdpbep86},
      {"dsd-pbep86", dsdpbep86},
      {"dsdsvwn", dsdsvwn},
      {"dsd-svwn", dsdsvwn},
      {"dodblyp", dodblyp},
      {"dod-blyp", dodblyp},
      {"dodpbe", dodpbe},
      {"dod-pbe", dodpbe},
      {"dodpbeb95", dodpbeb95},
      {"dod-pbeb95", dodpbeb95},
      {"dodpbep86", dodpbep86},
      {"dod-pbep86", dodpbep86},
      {"dodsvwn", dodsvwn},
      {"dod-svwn", dodsvwn},
      {"pbe02", pbe0_2},
      {"pbe0-2", pbe0_2},
      {"pbe0dh", pbe0_dh},
      {"pbe0-dh", pbe0_dh},
      {"hf-3c", hf3c},
      {"hf3c", hf3c},
      {"hf-3cv", hf3cv},
      {"hf3cv", hf3cv},
      {"pbeh3c", pbeh3c},
      {"pbeh-3c", pbeh3c},
      {"b973c", b973c},
      {"b97-3c", b973c},
      {"hsesol", hsesol},
      {"pwgga", pwgga},
      {"dftb3", dftb3},
      {"hcth120", hcth120},
      {"ptpss", ptpss},
      {"lc-wpbe", lcwpbe},
      {"lcwpbe", lcwpbe},
      {"bmk", bmk},
      {"b1b95", b1b95},
      {"pwb6k", pwb6k},
      {"otpss", otpss},
      {"ssb", ssb},
      {"revssb", revssb},
      {"pbesol", pbesol},
      {"hse06", hse06},
      {"pbexalpha", pbexalpha},
      {"pbehpbe", pbehpbe},
      {"hcth407", hcth407},
      {"n12", n12},
      {"pkzb", pkzb},
      {"thcth", thcth},
      {"tauhctc", thcth},
      {"m11l", m11l},
      {"mn15l", mn15l},
      {"mpwb1k", mpwb1k},
      {"mpw1kcis", mpw1kcis},
      {"mpwkcis1k", mpwkcis1k},
      {"pbeh1pbe", pbeh1pbe},
      {"pbe1kcis", pbe1kcis},
      {"b97-1", b97_1},
      {"b97-2", b97_2},
      {"b98", b98},
      {"hiss", hiss},
      {"hse03", hse03},
      {"revtpssh", revtpssh},
      {"tpss1kcis", tpss1kcis},
      {"m05", m05},
      {"m052x", m052x},
      {"m05-2x", m052x},
      {"m08hx", m08hx},
      {"m08-hx", m08hx},
      {"lcwhpbe", lcwhpbe},
      {"lc-whpbe", lcwhpbe},
      {"mn12l", mn12l},
      {"tauhcthhyb", tauhcthhyb},
      {"sogga11x", sogga11x},
      {"n12sx", n12sx},
      {"mn12sx", mn12sx},
      {"mn15", mn15},
      {"glyp", glyp},
      {"g-lyp", glyp},
      {"revpbe0dh", revpbe0dh},
      {"revpbe0-dh", revpbe0dh},
      {"revtpss0", revtpss0},
      {"revdodpbep86", revdodpbep86},
      {"revdod-pbep86", revdodpbep86},
      {"revdsdpbep86", revdsdpbep86},
      {"revdsd-pbep86", revdsdpbep86},
      {"revdsdpbe", revdsdpbe},
      {"revdsd-pbe", revdsdpbe},
      {"revdsdpbepbe", revdsdpbe},
      {"revdsd-pbepbe", revdsdpbe},
      {"revdsdblyp", revdsdblyp},
      {"revdsd-blyp", revdsdblyp},
      {"am05", am05},
      {"hse12", hse12},
      {"hse12s", hse12s},
      {"rscan", rscan},
      {"r2scan", r2scan},
      {"r2scanh", r2scanh},
      {"r2scan0", r2scan0},
      {"r2scan50", r2scan50},
      {"dftb3", dftb3},
      {"dftb_3ob", dftb3},
      {"dftb(3ob)", dftb3},
      {"dftb_mio", dftb_mio},
      {"dftb(mio)", dftb_mio},
      {"dftb_ob2", dftb_ob2},
      {"lc-dftb", dftb_ob2},
      {"dftb(ob2)", dftb_ob2},
      {"dftb_matsci", dftb_matsci},
      {"dftb(matsci)", dftb_matsci},
      {"dftb_pbc", dftb_pbc},
      {"dftb(pbc)", dftb_pbc},
  };

  std::string func = name;
  transform(name.begin(), name.end(), func.begin(), ::tolower);
  auto iter = dfuncs.find(func);
  if (iter != dfuncs.end()) {
    return iter->second;
  }
  return none;
};

dparam get_d4eeqbjatm_2019_parameter(dfunc num) {
  dparam par;
  double s6{1.0}, s8{0.0}, s10{0.0}, s9{1.0}, a1{0.0}, a2{0.0};
  int alp{16};
  switch (num) {
    default:
      break;
    case am05: // (SAW211021)
      s8 = 1.71885838;
      a1 = 0.47901431;
      a2 = 5.96771581;
      break;
      // Fitset: MD= -0.28899 MAD= 0.52215 RMSD= 0.93584
    case b1b95:
      s8 = 1.27701162;
      a1 = 0.40554715;
      a2 = 4.63323074;
      break;
      // Fitset: MD= 0.22852 MAD= 0.35189 RMSD= 0.46982
    case b1lyp:
      s8 = 1.98553711;
      a1 = 0.39309040;
      a2 = 4.55465145;
      break;
      // Fitset: MD= -0.04797 MAD= 0.25597 RMSD= 0.38778
    case b1p:
      s8 = 3.36115015;
      a1 = 0.48665293;
      a2 = 5.05219572;
      break;
      // Fitset: MD= -0.01406 MAD= 0.27441 RMSD= 0.47328
    case b1pw:
      s8 = 3.02227550;
      a1 = 0.47396846;
      a2 = 4.49845309;
      break;
      // Fitset: MD= 0.10485 MAD= 0.32175 RMSD= 0.48508
    case b2gpplyp:
      s6 = 0.5600;
      s8 = 0.94633372;
      a1 = 0.42907301;
      a2 = 5.18802602;
      break;
      // Fitset: MD= -0.05248 MAD= 0.18110 RMSD= 0.27365
    case b2plyp:
      s6 = 0.6400;
      s8 = 1.16888646;
      a1 = 0.44154604;
      a2 = 4.73114642;
      break;
      // Fitset: MD= -0.03761 MAD= 0.18247 RMSD= 0.27109
    case b3lyp:
      s8 = 2.02929367;
      a1 = 0.40868035;
      a2 = 4.53807137;
      break;
      // Fitset: MD= -0.05892 MAD= 0.26117 RMSD= 0.40531
    case b3p:
      s8 = 3.08822155;
      a1 = 0.47324238;
      a2 = 4.98682134;
      break;
      // Fitset: MD= -0.02970 MAD= 0.26962 RMSD= 0.46761
    case b3pw:
      s8 = 2.88364295;
      a1 = 0.46990860;
      a2 = 4.51641422;
      break;
      // Fitset: MD= 0.06643 MAD= 0.29151 RMSD= 0.45541
    case b97:
      s8 = 0.87854260;
      a1 = 0.29319126;
      a2 = 4.51647719;
      break;
      // Fitset: MD= -0.13017 MAD= 0.24778 RMSD= 0.36116
    case bhlyp:
      s8 = 1.65281646;
      a1 = 0.27263660;
      a2 = 5.48634586;
      break;
      // Fitset: MD= -0.15832 MAD= 0.34132 RMSD= 0.57342
    case blyp:
      s8 = 2.34076671;
      a1 = 0.44488865;
      a2 = 4.09330090;
      break;
      // Fitset: MD= 0.04801 MAD= 0.28161 RMSD= 0.38321
    case bpbe:
      s8 = 3.64405246;
      a1 = 0.52905620;
      a2 = 4.11311891;
      break;
      // Fitset: MD= 0.19316 MAD= 0.41912 RMSD= 0.60452
    case bp:
      s8 = 3.35497927;
      a1 = 0.43645861;
      a2 = 4.92406854;
      break;
      // Fitset: MD= 0.08252 MAD= 0.32681 RMSD= 0.47063
    case bpw:
      s8 = 3.24571506;
      a1 = 0.50050454;
      a2 = 4.12346483;
      break;
      // Fitset: MD= 0.20607 MAD= 0.41941 RMSD= 0.59589
    case camb3lyp:
      s8 = 1.66041301;
      a1 = 0.40267156;
      a2 = 5.17432195;
      break;
      // Fitset: MD= -0.19675 MAD= 0.34901 RMSD= 0.59087
   case dftb3: // (SAW191202)
      s6 = 1.0;
      s8 = 0.6635015;
      a1 = 0.5523240;
      a2 = 4.3537076;
      break;
   case dftb_matsci: // (SAW191202)
      s6 = 1.0;
      s8 = 3.3157614;
      a1 = 0.4826330;
      a2 = 5.3811976;
      break;
   case dftb_mio: // (SAW191202)
      s6 = 1.0;
      s8 = 1.2916225;
      a1 = 0.5965326;
      a2 = 4.8778602;
      break;
   case dftb_ob2: // (SAW191202)
      s6 = 1.0;
      s8 = 2.9692689;
      a1 = 0.6068916;
      a2 = 5.4476789;
      break;
   case dftb_pbc: // (SAW191202)
      s6 = 1.0;
      s8 = 2.1667394;
      a1 = 0.5646391;
      a2 = 4.9576353;
      break;
    case dodblyp:
      s6 = 0.4700;
      s8 = 1.31146043;
      a1 = 0.43407294;
      a2 = 4.27914360;
      break;
      // Fitset: MD= 0.03323 MAD= 0.13858 RMSD= 0.20861
    case dodpbeb95:
      s6 = 0.5600;
      s8 = 0.01574635;
      a1 = 0.43745720;
      a2 = 3.69180763;
      break;
      // Fitset: MD= 0.03704 MAD= 0.13343 RMSD= 0.18278
    case dodpbe:
      s6 = 0.4800;
      s8 = 0.92051454;
      a1 = 0.43037052;
      a2 = 4.38067238;
      break;
      // Fitset: MD= 0.01065 MAD= 0.13414 RMSD= 0.21424
    case dodpbep86:
      s6 = 0.4600;
      s8 = 0.71405681;
      a1 = 0.42408665;
      a2 = 4.52884439;
      break;
      // Fitset: MD= -0.03740 MAD= 0.12467 RMSD= 0.18127
    case dodsvwn:
      s6 = 0.4200;
      s8 = 0.94500207;
      a1 = 0.47449026;
      a2 = 5.05316093;
      break;
      // Fitset: MD= -0.07427 MAD= 0.16970 RMSD= 0.25286
    case dsdblyp:
      s6 = 0.5400;
      s8 = 0.63018237;
      a1 = 0.47591835;
      a2 = 4.73713781;
      break;
      // Fitset: MD= -0.01981 MAD= 0.14823 RMSD= 0.21530
    case dsdpbeb95:
      s6 = 0.5400;
      s8 = -0.14668670;
      a1 = 0.46394587;
      a2 = 3.64913860;
      break;
      // Fitset: MD= 0.02996 MAD= 0.12414 RMSD= 0.16860
    case dsdpbe:
      s6 = 0.4500;
      s8 = 0.70584116;
      a1 = 0.45787085;
      a2 = 4.44566742;
      break;
      // Fitset: MD= 0.00866 MAD= 0.13406 RMSD= 0.21380
    case dsdpbep86:
      s6 = 0.4700;
      s8 = 0.37586675;
      a1 = 0.53698768;
      a2 = 5.13022435;
      break;
      // Fitset: MD= -0.05273 MAD= 0.14259 RMSD= 0.21271
    case dsdsvwn:
      s6 = 0.4100;
      s8 = 0.72914436;
      a1 = 0.51347412;
      a2 = 5.11858541;
      break;
      // Fitset: MD= -0.08974 MAD= 0.32285 RMSD= 0.43146
    case glyp:
      s8 = 4.23798924;
      a1 = 0.38426465;
      a2 = 4.38412863;
      break;
      // Fitset: MD= 0.63466 MAD= 0.89568 RMSD= 1.11309
    case hf:
      s8 = 1.61679827;
      a1 = 0.44959224;
      a2 = 3.35743605;
      break;
      // Fitset: MD= -0.02597 MAD= 0.34732 RMSD= 0.49719
    case(hse03): // (SAW211107)
      s6 = 1.0;
      s8 = 1.19812280;
      a1 = 0.38662939;
      a2 = 5.22925796;
      break;
    case(hse06): // (SAW211107)
      s6 = 1.0;
      s8 = 1.19528249;
      a1 = 0.38663183;
      a2 = 5.19133469;
      break;
    case(hse12): // (SAW211107)
      s6 = 1.0;
      s8 = 1.23500792;
      a1 = 0.39226921;
      a2 = 5.22036266;
      break;
    case(hse12s): // (SAW211107)
      s6 = 1.0;
      s8 = 1.23767762;
      a1 = 0.39989137;
      a2 = 5.34809245;
      break;
    case(hsesol): // (SAW211107)
      s6 = 1.0;
      s8 = 1.82207807;
      a1 = 0.45646268;
      a2 = 5.59662251;
      break;
    case lb94:
      s8 = 2.59538499;
      a1 = 0.42088944;
      a2 = 3.28193223;
      break;
      // Fitset: MD= 0.31701 MAD= 0.53196 RMSD= 0.74553
    case lcblyp:
      s8 = 1.60344180;
      a1 = 0.45769839;
      a2 = 7.86924893;
      break;
      // Fitset: MD= -0.39724 MAD= 0.72327 RMSD= 1.18218
    case lh07ssvwn:
      s8 = 3.16675531;
      a1 = 0.35965552;
      a2 = 4.31947614;
      break;
      // Fitset: MD= 0.32224 MAD= 0.59006 RMSD= 0.86272
    case lh07tsvwn:
      s8 = 2.09333001;
      a1 = 0.35025189;
      a2 = 4.34166515;
      break;
      // Fitset: MD= 0.24243 MAD= 0.43497 RMSD= 0.61671
    case lh12ctssifpw92:
      s8 = 2.68467610;
      a1 = 0.34190416;
      a2 = 3.91039666;
      break;
      // Fitset: MD= 0.55106 MAD= 0.80783 RMSD= 1.11048
    case lh12ctssirpw92:
      s8 = 2.48973402;
      a1 = 0.34026075;
      a2 = 3.96948081;
      break;
      // Fitset: MD= 0.47785 MAD= 0.71188 RMSD= 0.98422
    case lh14tcalpbe:
      s8 = 1.28130770;
      a1 = 0.38822021;
      a2 = 4.92501211;
      break;
      // Fitset: MD= -0.02105 MAD= 0.22968 RMSD= 0.36045
    case lh20t: // (10.1021/acs.jctc.0c00498)
      s6 = 1.000;
      s8 = 0.113;
      a1 = 0.479;
      a2 = 4.635;
      break;
    case m06:
      s8 = 0.16366729;
      a1 = 0.53456413;
      a2 = 6.06192174;
      break;
      // Fitset: MD= 0.01788 MAD= 0.24914 RMSD= 0.38604
    case m06l:
      s8 = 0.59493760;
      a1 = 0.71422359;
      a2 = 6.35314182;
      break;
      // Fitset: MD= 0.08395 MAD= 0.24888 RMSD= 0.34879
    case(mn12sx): // (SAW211021)
      s8 = 0.85964873;
      a1 = 0.62662681;
      a2 = 5.62088906;
      break;
      // Fitset: MD= 0.16131 MAD= 0.34142 RMSD= 0.47113
    case mpw1b95:
      s8 = 0.50093024;
      a1 = 0.41585097;
      a2 = 4.99154869;
      break;
      // Fitset: MD= 0.00585 MAD= 0.15695 RMSD= 0.21297
    case mpw1lyp:
      s8 = 1.15591153;
      a1 = 0.25603493;
      a2 = 5.32083895;
      break;
      // Fitset: MD= -0.26979 MAD= 0.41542 RMSD= 0.60678
    case mpw1pw:
      s8 = 1.80841716;
      a1 = 0.42961819;
      a2 = 4.68892341;
      break;
      // Fitset: MD= -0.08840 MAD= 0.26815 RMSD= 0.45231
    case mpw2plyp:
      s6 = 0.7500;
      s8 = 0.45788846;
      a1 = 0.42997704;
      a2 = 5.07650682;
      break;
      // Fitset: MD= -0.18921 MAD= 0.30115 RMSD= 0.44049
    case mpwb1k:
      s8 = 0.57338313;
      a1 = 0.44687975;
      a2 = 5.21266777;
      break;
      // Fitset: MD= -0.00870 MAD= 0.17226 RMSD= 0.23614
    case mpwlyp:
      s8 = 1.25842942;
      a1 = 0.25773894;
      a2 = 5.02319542;
      break;
      // Fitset: MD= -0.24426 MAD= 0.39145 RMSD= 0.54503
    case mpwpw:
      s8 = 1.82596836;
      a1 = 0.34526745;
      a2 = 4.84620734;
      break;
      // Fitset: MD= -0.06278 MAD= 0.27913 RMSD= 0.43988
    case o3lyp:
      s8 = 1.75762508;
      a1 = 0.10348980;
      a2 = 6.16233282;
      break;
      // Fitset: MD= -0.19268 MAD= 0.38577 RMSD= 0.62168
    case olyp:
      s8 = 2.74836820;
      a1 = 0.60184498;
      a2 = 2.53292167;
      break;
      // Fitset: MD= 0.12352 MAD= 0.37113 RMSD= 0.58291
    case opbe:
      s8 = 3.06917417;
      a1 = 0.68267534;
      a2 = 2.22849018;
      break;
      // Fitset: MD= 0.26699 MAD= 0.55308 RMSD= 0.85023
    case pbe0_2:
      s6 = 0.5000;
      s8 = 0.64299082;
      a1 = 0.76542115;
      a2 = 5.78578675;
      break;
      // Fitset: MD= -0.04260 MAD= 0.21186 RMSD= 0.34045
    case pbe0:
      s8 = 1.20065498;
      a1 = 0.40085597;
      a2 = 5.02928789;
      break;
      // Fitset: MD= -0.17892 MAD= 0.30557 RMSD= 0.51050
    case pbe0_dh:
      s6 = 0.8750;
      s8 = 0.96811578;
      a1 = 0.47592488;
      a2 = 5.08622873;
      break;
      // Fitset: MD= -0.13857 MAD= 0.27919 RMSD= 0.47256
    case pbe:
      s8 = 0.95948085;
      a1 = 0.38574991;
      a2 = 4.80688534;
      break;
      // Fitset: MD= -0.20544 MAD= 0.33635 RMSD= 0.51168
    case(pbesol): // (SAW211021)
      s8 = 1.71885698;
      a1 = 0.47901421;
      a2 = 5.96771589;
      break;
      // Fitset: MD= -0.28899 MAD= 0.52215 RMSD= 0.93584
    case pw1pw:
      s8 = 0.96850170;
      a1 = 0.42427511;
      a2 = 5.02060636;
      break;
      // Fitset: MD= -0.27325 MAD= 0.42206 RMSD= 0.64119
    case pw6b95:
      s8 = -0.31926054;
      a1 = 0.04142919;
      a2 = 5.84655608;
      break;
      // Fitset: MD= -0.04767 MAD= 0.14330 RMSD= 0.18958
    case pw86pbe:
      s8 = 1.21362856;
      a1 = 0.40510366;
      a2 = 4.66737724;
      break;
      // Fitset: MD= -0.11505 MAD= 0.24691 RMSD= 0.38101
    case pw91:
      s8 = 0.77283111;
      a1 = 0.39581542;
      a2 = 4.93405761;
      break;
      // Fitset: MD= -0.33019 MAD= 0.48611 RMSD= 0.68110
    case pwp1:
      s8 = 0.60492565;
      a1 = 0.46855837;
      a2 = 5.76921413;
      break;
      // Fitset: MD= -0.35321 MAD= 0.54026 RMSD= 0.86629
    case pwpb95:
      s6 = 0.8200;
      s8 = -0.34639127;
      a1 = 0.41080636;
      a2 = 3.83878274;
      break;
      // Fitset: MD= 0.02143 MAD= 0.13040 RMSD= 0.17599
    case pwp:
      s8 = 0.32801227;
      a1 = 0.35874687;
      a2 = 6.05861168;
      break;
      // Fitset: MD= -0.42482 MAD= 0.62607 RMSD= 0.91840
    case(revdsdpbep86):
      s6 = 0.5132;
      s8 = 0.00000000;
      a1 = 0.44000000;
      a2 = 3.60000000;
      break;
    case(revdsdpbe):
      s6 = 0.6706;
      s8 = 0.00000000;
      a1 = 0.40000000;
      a2 = 3.60000000;
      break;
    case(revdsdblyp):
      s6 = 0.6141;
      s8 = 0.00000000;
      a1 = 0.38000000;
      a2 = 3.52000000;
      break;
    case(revdodpbep86):
      s6 = 0.5552;
      s8 = 0.00000000;
      a1 = 0.44000000;
      a2 = 3.60000000;
      break;
    case revpbe0:
      s8 = 1.57185414;
      a1 = 0.38705966;
      a2 = 4.11028876;
      break;
      // Fitset: MD= 0.02724 MAD= 0.21587 RMSD= 0.36040
    case revpbe0dh:
      s6 = 0.8750;
      s8 = 1.24456037;
      a1 = 0.36730560;
      a2 = 4.71126482;
      break;
      // Fitset: MD= -0.01089 MAD= 0.20910 RMSD= 0.33564
    case revpbe38:
      s8 = 1.66597472;
      a1 = 0.39476833;
      a2 = 4.39026628;
      break;
      // Fitset: MD= -0.01326 MAD= 0.22598 RMSD= 0.36210
    case revpbe:
      s8 = 1.74676530;
      a1 = 0.53634900;
      a2 = 3.07261485;
      break;
      // Fitset: MD= 0.05649 MAD= 0.25212 RMSD= 0.40863
    case revtpss0:
      s8 = 1.54664499;
      a1 = 0.45890964;
      a2 = 4.78426405;
      break;
      // Fitset: MD= -0.05298 MAD= 0.19965 RMSD= 0.32081
    case revtpss:
      s8 = 1.53089454;
      a1 = 0.44880597;
      a2 = 4.64042317;
      break;
      // Fitset: MD= -0.01904 MAD= 0.19568 RMSD= 0.29618
    case revtpssh:
      s8 = 1.52740307;
      a1 = 0.45161957;
      a2 = 4.70779483;
      break;
      // Fitset: MD= -0.03731 MAD= 0.19133 RMSD= 0.29091
    case rpbe:
      s8 = 1.31183787;
      a1 = 0.46169493;
      a2 = 3.15711757;
      break;
      // Fitset: MD= -0.07156 MAD= 0.26348 RMSD= 0.38671
    case rpw86pbe:
      s8 = 1.12624034;
      a1 = 0.38151218;
      a2 = 4.75480472;
      break;
      // Fitset: MD= -0.12740 MAD= 0.26294 RMSD= 0.40614
    case scan:
      s8 = 1.46126056;
      a1 = 0.62930855;
      a2 = 6.31284039;
      break;
    case(rscan): // (10.1063/5.0041008)
      s8 = 0.87728975;
      a1 = 0.49116966;
      a2 = 5.75859346;
      break;
   case(r2scan): // (10.1063/5.0041008)
      s8 = 0.60187490;
      a1 = 0.51559235;
      a2 = 5.77342911;
      break;
   case(r2scanh):
      s6 = 1.0;
      s8 = 0.8324;
      a1 = 0.4944;
      a2 = 5.9019;
      break;
   case(r2scan0):
      s6 = 1.0;
      s8 = 0.8992;
      a1 = 0.4778;
      a2 = 5.8779;
      break;
   case(r2scan50):
      s6 = 1.0;
      s8 = 1.0471;
      a1 = 0.4574;
      a2 = 5.8969;
      break;
      // Fitset: MD= -0.13170 MAD= 0.28640 RMSD= 0.51183
    case tpss0:
      s8 = 1.62438102;
      a1 = 0.40329022;
      a2 = 4.80537871;
      break;
      // Fitset: MD= -0.09569 MAD= 0.26733 RMSD= 0.44767
    case tpss:
      s8 = 1.76596355;
      a1 = 0.42822303;
      a2 = 4.54257102;
      break;
      // Fitset: MD= -0.09296 MAD= 0.27505 RMSD= 0.42537
    case tpssh:
      s8 = 1.85897750;
      a1 = 0.44286966;
      a2 = 4.60230534;
      break;
      // Fitset: MD=  0.02238 MAD= 0.16042 RMSD= 0.33519
    case b97d: // (SAW201029)
      s8 = 1.69460052;
      a1 = 0.28904684;
      a2 = 4.13407323;
      break;
      // Fitset: MD= -0.09858 MAD= 0.26757 RMSD= 0.42380
    case wb97: // (SAW190103)
      s8 = 6.55792598;
      a1 = 0.76666802;
      a2 = 8.36027334;
      break;
      // Fitset: MD= -0.12779 MAD= 0.36152 RMSD= 0.49991
    case wb97x: // (SAW190103)
      s8 = -0.07519516;
      a1 = 0.45094893;
      a2 = 6.78425255;
      break;
      // S22x5: MD= 0.05 MAD= 0.16 RMSD= 0.22
      // S66x8: MD= 0.06 MAD= 0.16 RMSD= 0.21
      // NCI10: MD= 0.08 MAD= 0.15 RMSD= 0.25
    case b97m: // (10.1002/jcc.26411)
      s8 = 0.6633;
      a1 = 0.4288;
      a2 = 3.9935;
      break;
      // S22x5: MD= 0.03 MAD= 0.12 RMSD= 0.18
      // S66x8: MD= 0.09 MAD= 0.17 RMSD= 0.22
      // NCI10: MD= 0.09 MAD= 0.15 RMSD= 0.32
    case wb97m: // (10.1002/jcc.26411)
      s8 = 0.7761;
      a1 = 0.7514;
      a2 = 2.7099;
      break;
      // Fitset: MD= -0.20216 MAD= 0.34696 RMSD= 0.53641
    case x3lyp:
      s8 = 1.54701429;
      a1 = 0.20318443;
      a2 = 5.61852648;
      break;
      // Fitset: MD= -0.15607 MAD= 0.31342 RMSD= 0.49546
    case xlyp:
      s8 = 1.62972054;
      a1 = 0.11268673;
      a2 = 5.40786417;
      break;
      // Fitset: MD= -0.03900 MAD= 0.27562 RMSD= 0.38491
  }

  par.s6 = s6;
  par.s8 = s8;
  par.s9 = s9;
  par.s10 = s10;
  par.a1 = a1;
  par.a2 = a2;
  par.alp = alp;

  return par;
}

dparam get_d4eeqbjmbd_2019_parameter(dfunc num) {
  dparam par;
  double s6{1.0}, s8{0.0}, s10{0.0}, s9{1.0}, a1{0.0}, a2{0.0};
  int alp{16};
  switch (num) {
    default:
      break;
    case b1b95:
      s8 = 1.19549420;
      a1 = 0.39241474;
      a2 = 4.60397611;
      break;
      // Fitset: MD= 0.21329 MAD= 0.33289 RMSD= 0.44693
    case b1lyp:
      s8 = 1.94609514;
      a1 = 0.38643351;
      a2 = 4.54135968;
      break;
      // Fitset: MD= -0.06493 MAD= 0.25607 RMSD= 0.39776
    case b1p:
      s8 = 3.38693011;
      a1 = 0.48478615;
      a2 = 5.04361224;
      break;
      // Fitset: MD= -0.02348 MAD= 0.27543 RMSD= 0.48014
    case b1pw:
      s8 = 2.98402204;
      a1 = 0.46862950;
      a2 = 4.48637849;
      break;
      // Fitset: MD= 0.09181 MAD= 0.31824 RMSD= 0.48100
    case b2gpplyp:
      s8 = 1.00494214;
      a1 = 0.42447353;
      a2 = 5.19461329;
      break;
      // Fitset: MD= -0.06339 MAD= 0.18624 RMSD= 0.28316
    case b2plyp:
      s6 = 0.6400;
      s8 = 1.15117773;
      a1 = 0.42666167;
      a2 = 4.73635790;
      break;
      // Fitset: MD= -0.05031 MAD= 0.18506 RMSD= 0.28010
    case b3lyp:
      s8 = 2.00246246;
      a1 = 0.40276191;
      a2 = 4.52778320;
      break;
      // Fitset: MD= -0.07554 MAD= 0.26205 RMSD= 0.41586
    case b3p:
      s8 = 3.14456298;
      a1 = 0.47187947;
      a2 = 4.98624258;
      break;
      // Fitset: MD= -0.04085 MAD= 0.27080 RMSD= 0.47542
    case b3pw:
      s8 = 2.85656268;
      a1 = 0.46491801;
      a2 = 4.50601452;
      break;
      // Fitset: MD= 0.05302 MAD= 0.28885 RMSD= 0.45409
    case b97:
      s8 = 0.81171211;
      a1 = 0.28461283;
      a2 = 4.48691468;
      break;
      // Fitset: MD= -0.15506 MAD= 0.26269 RMSD= 0.38008
    case bhlyp:
      s8 = 1.68082973;
      a1 = 0.26835837;
      a2 = 5.48847218;
      break;
      // Fitset: MD= -0.17297 MAD= 0.34819 RMSD= 0.58920
    case blyp:
      s8 = 2.33971306;
      a1 = 0.44733688;
      a2 = 4.06583931;
      break;
      // Fitset: MD= 0.02464 MAD= 0.27084 RMSD= 0.37647
    case bpbe:
      s8 = 3.65322996;
      a1 = 0.49933501;
      a2 = 4.24294852;
      break;
      // Fitset: MD= 0.18169 MAD= 0.41723 RMSD= 0.60190
    case bp:
      s8 = 3.33728176;
      a1 = 0.43220330;
      a2 = 4.91443061;
      break;
      // Fitset: MD= 0.07108 MAD= 0.32331 RMSD= 0.46959
    case bpw:
      s8 = 3.23137432;
      a1 = 0.49955226;
      a2 = 4.10411084;
      break;
      // Fitset: MD= 0.18835 MAD= 0.41014 RMSD= 0.58149
    case camb3lyp:
      s8 = 1.74407961;
      a1 = 0.40137870;
      a2 = 5.18731225;
      break;
      // Fitset: MD= -0.20897 MAD= 0.35917 RMSD= 0.60661
    case dodblyp:
      s6 = 0.4700;
      s8 = 1.17809956;
      a1 = 0.40252428;
      a2 = 4.25096555;
      break;
      // Fitset: MD= 0.01949 MAD= 0.13138 RMSD= 0.19558
    case dodpbeb95:
      s6 = 0.5400;
      s8 = -0.15702803;
      a1 = 0.30629389;
      a2 = 3.69170956;
      break;
      // Fitset: MD= -0.00076 MAD= 0.09552 RMSD= 0.12939
    case dodpbe:
      s6 = 0.4800;
      s8 = 0.83908332;
      a1 = 0.40655901;
      a2 = 4.33601239;
      break;
      // Fitset: MD= -0.00580 MAD= 0.12659 RMSD= 0.20494
    case dodpbep86:
      s6 = 0.4600;
      s8 = 0.68309910;
      a1 = 0.40600975;
      a2 = 4.50011772;
      break;
      // Fitset: MD= -0.05518 MAD= 0.12509 RMSD= 0.18597
    case dodsvwn:
      s6 = 0.4200;
      s8 = 1.01890345;
      a1 = 0.46167459;
      a2 = 5.11121382;
      break;
      // Fitset: MD= -0.08458 MAD= 0.17447 RMSD= 0.26263
    case dsdblyp:
      s6 = 0.5400;
      s8 = 0.65438817;
      a1 = 0.46549574;
      a2 = 4.73449899;
      break;
      // Fitset: MD= -0.03244 MAD= 0.15073 RMSD= 0.22004
    case dsdpbeb95:
      s6 = 0.5400;
      s8 = -0.24336862;
      a1 = 0.32697409;
      a2 = 3.69767540;
      break;
      // Fitset: MD= -0.00830 MAD= 0.09184 RMSD= 0.12546
    case dsdpbe:
      s6 = 0.4500;
      s8 = 0.66116783;
      a1 = 0.43565915;
      a2 = 4.41110670;
      break;
      // Fitset: MD= -0.00698 MAD= 0.12881 RMSD= 0.20652
    case dsdpbep86:
      s6 = 0.4700;
      s8 = 0.51157821;
      a1 = 0.53889789;
      a2 = 5.18645943;
      break;
      // Fitset: MD= -0.06077 MAD= 0.14727 RMSD= 0.21987
    case dsdsvwn:
      s6 = 0.4100;
      s8 = 0.90084457;
      a1 = 0.51106529;
      a2 = 5.22490148;
      break;
      // Fitset: MD= -0.08917 MAD= 0.32295 RMSD= 0.43170
    case glyp:
      s8 = 3.83861584;
      a1 = 0.36343954;
      a2 = 4.32875183;
      break;
      // Fitset: MD= 0.63864 MAD= 0.89437 RMSD= 1.11594
    case hf:
      s8 = 1.46001146;
      a1 = 0.43186901;
      a2 = 3.34116014;
      break;
      // Fitset: MD= -0.05320 MAD= 0.34504 RMSD= 0.50301
    case lb94:
      s8 = 2.36461524;
      a1 = 0.41518379;
      a2 = 3.19365471;
      break;
      // Fitset: MD= 0.28794 MAD= 0.49917 RMSD= 0.70755
    case lcblyp:
      s8 = 2.40109962;
      a1 = 0.47867438;
      a2 = 8.01038424;
      break;
      // Fitset: MD= -0.39664 MAD= 0.72492 RMSD= 1.18496
    case lh07ssvwn:
      s8 = 2.92498406;
      a1 = 0.34173988;
      a2 = 4.28404951;
      break;
      // Fitset: MD= 0.31275 MAD= 0.58894 RMSD= 0.85534
    case lh07tsvwn:
      s8 = 1.95389300;
      a1 = 0.33511515;
      a2 = 4.31853958;
      break;
      // Fitset: MD= 0.22401 MAD= 0.42077 RMSD= 0.59849
    case lh12ctssifpw92:
      s8 = 2.41356607;
      a1 = 0.31391316;
      a2 = 3.88935769;
      break;
      // Fitset: MD= 0.53322 MAD= 0.78975 RMSD= 1.08954
    case lh12ctssirpw92:
      s8 = 2.24917162;
      a1 = 0.31446575;
      a2 = 3.95070925;
      break;
      // Fitset: MD= 0.45858 MAD= 0.69302 RMSD= 0.96196
    case lh14tcalpbe:
      s8 = 1.27677253;
      a1 = 0.38128670;
      a2 = 4.91698883;
      break;
      // Fitset: MD= -0.03475 MAD= 0.22645 RMSD= 0.36554
    case m06:
      s8 = 0.22948274;
      a1 = 0.52927285;
      a2 = 6.06516782;
      break;
      // Fitset: MD= 0.01376 MAD= 0.24790 RMSD= 0.38566
    case m06l:
      s8 = 0.40077779;
      a1 = 0.69611405;
      a2 = 6.29092087;
      break;
      // Fitset: MD= 0.08204 MAD= 0.24719 RMSD= 0.34728
    case mpw1b95:
      s8 = 0.53791835;
      a1 = 0.41016913;
      a2 = 4.99284176;
      break;
      // Fitset: MD= -0.00740 MAD= 0.15464 RMSD= 0.20903
    case mpw1lyp:
      s8 = 1.19986100;
      a1 = 0.25502469;
      a2 = 5.32301304;
      break;
      // Fitset: MD= -0.28762 MAD= 0.43446 RMSD= 0.62777
    case mpw1pw:
      s8 = 1.80656973;
      a1 = 0.42456967;
      a2 = 4.68132317;
      break;
      // Fitset: MD= -0.10273 MAD= 0.27455 RMSD= 0.46383
    case mpw2plyp:
      s6 = 0.7500;
      s8 = 0.61161179;
      a1 = 0.43748316;
      a2 = 5.12540364;
      break;
      // Fitset: MD= -0.20058 MAD= 0.31617 RMSD= 0.45768
    case mpwb1k:
      s8 = 0.62221146;
      a1 = 0.44216745;
      a2 = 5.21324659;
      break;
      // Fitset: MD= -0.01872 MAD= 0.17358 RMSD= 0.23758
    case mpwlyp:
      s8 = 1.18243337;
      a1 = 0.38968985;
      a2 = 4.30835285;
      break;
      // Fitset: MD= -0.27161 MAD= 0.42078 RMSD= 0.58498
    case mpwpw:
      s8 = 1.79674014;
      a1 = 0.33870479;
      a2 = 4.83442213;
      break;
      // Fitset: MD= -0.07933 MAD= 0.28408 RMSD= 0.44655
    case o3lyp:
      s8 = 1.77793802;
      a1 = 0.09961745;
      a2 = 6.16089304;
      break;
      // Fitset: MD= -0.20876 MAD= 0.39753 RMSD= 0.63678
    case olyp:
      s8 = 2.58717041;
      a1 = 0.59759271;
      a2 = 2.48760353;
      break;
      // Fitset: MD= 0.09561 MAD= 0.37723 RMSD= 0.57321
    case opbe:
      s8 = 2.93544102;
      a1 = 0.67903933;
      a2 = 2.19810071;
      break;
      // Fitset: MD= 0.23963 MAD= 0.55856 RMSD= 0.83684
    case pbe0_2:
      s6 = 0.5000;
      s8 = 0.98834859;
      a1 = 0.77911062;
      a2 = 5.90389569;
      break;
      // Fitset: MD= -0.04530 MAD= 0.21377 RMSD= 0.34234
    case pbe0:
      s8 = 1.26829475;
      a1 = 0.39907098;
      a2 = 5.03951304;
      break;
      // Fitset: MD= -0.19236 MAD= 0.31907 RMSD= 0.52730
    case pbe0_dh:
      s6 = 0.8750;
      s8 = 1.19306002;
      a1 = 0.46106784;
      a2 = 5.25210480;
      break;
      // Fitset: MD= -0.15041 MAD= 0.29062 RMSD= 0.48405
    case pbe:
      s8 = 0.99924614;
      a1 = 0.38142528;
      a2 = 4.81839284;
      break;
      // Fitset: MD= -0.22162 MAD= 0.35068 RMSD= 0.52976
    case pw1pw:
      s8 = 1.09759050;
      a1 = 0.42759830;
      a2 = 5.04559572;
      break;
      // Fitset: MD= -0.28594 MAD= 0.43822 RMSD= 0.65926
    case pw6b95:
      s8 = -0.31629935;
      a1 = 0.03999357;
      a2 = 5.83690254;
      break;
      // Fitset: MD= -0.07192 MAD= 0.15341 RMSD= 0.20230
    case pw86pbe:
      s8 = 1.22842987;
      a1 = 0.39998824;
      a2 = 4.66739111;
      break;
      // Fitset: MD= -0.13164 MAD= 0.25352 RMSD= 0.39362
    case pw91:
      s8 = 0.81406882;
      a1 = 0.34094706;
      a2 = 5.18568823;
      break;
      // Fitset: MD= -0.34117 MAD= 0.49844 RMSD= 0.69615
    case pwp1:
      s8 = 0.95936222;
      a1 = 0.48552982;
      a2 = 5.84956411;
      break;
      // Fitset: MD= -0.35812 MAD= 0.55157 RMSD= 0.87759
    case pwpb95:
      s6 = 0.8200;
      s8 = -0.46453780;
      a1 = 0.29884136;
      a2 = 3.87641255;
      break;
      // Fitset: MD= -0.01702 MAD= 0.10667 RMSD= 0.14619
    case pwp:
      s8 = 0.66056055;
      a1 = 0.37768052;
      a2 = 6.14787138;
      break;
      // Fitset: MD= -0.42973 MAD= 0.63830 RMSD= 0.93118
    case revpbe0:
      s8 = 1.47198256;
      a1 = 0.37471756;
      a2 = 4.08904369;
      break;
      // Fitset: MD= 0.00549 MAD= 0.21492 RMSD= 0.35464
    case revpbe0dh:
      s6 = 0.8750;
      s8 = 1.22494188;
      a1 = 0.35904781;
      a2 = 4.70216012;
      break;
      // Fitset: MD= -0.02829 MAD= 0.21126 RMSD= 0.33843
    case revpbe38:
      s8 = 1.60423529;
      a1 = 0.38938475;
      a2 = 4.35557832;
      break;
      // Fitset: MD= -0.03169 MAD= 0.22643 RMSD= 0.36450
    case revpbe:
      s8 = 1.62543693;
      a1 = 0.54031831;
      a2 = 2.97965648;
      break;
      // Fitset: MD= 0.02784 MAD= 0.24691 RMSD= 0.39242
    case revtpss0:
      s8 = 1.55321888;
      a1 = 0.45355319;
      a2 = 4.77588598;
      break;
      // Fitset: MD= -0.06532 MAD= 0.20465 RMSD= 0.32888
    case revtpss:
      s8 = 1.51858035;
      a1 = 0.44243222;
      a2 = 4.62881620;
      break;
      // Fitset: MD= -0.03308 MAD= 0.19729 RMSD= 0.29587
    case revtpssh:
      s8 = 1.52542064;
      a1 = 0.44570207;
      a2 = 4.69883717;
      break;
      // Fitset: MD= -0.05069 MAD= 0.19530 RMSD= 0.29596
    case rpbe:
      s8 = 1.11793696;
      a1 = 0.44632488;
      a2 = 3.08890917;
      break;
      // Fitset: MD= -0.10098 MAD= 0.27012 RMSD= 0.39113
    case rpw86pbe:
      s8 = 1.13795871;
      a1 = 0.37636536;
      a2 = 4.75236384;
      break;
      // Fitset: MD= -0.14431 MAD= 0.26983 RMSD= 0.41992
    case scan:
      s8 = 1.75408315;
      a1 = 0.63571334;
      a2 = 6.35690748;
      break;
      // Fitset: MD= -0.13418 MAD= 0.28800 RMSD= 0.51482
    case tpss0:
      s8 = 1.66752698;
      a1 = 0.40074746;
      a2 = 4.80927196;
      break;
      // Fitset: MD= -0.11178 MAD= 0.27563 RMSD= 0.45972
    case tpss:
      s8 = 1.91130849;
      a1 = 0.43332851;
      a2 = 4.56986797;
      break;
      // Fitset: MD= -0.11454 MAD= 0.28494 RMSD= 0.43743
    case tpssh:
      s8 = 1.88783525;
      a1 = 0.43968167;
      a2 = 4.60342700;
      break;
      // Fitset: MD= -0.11391 MAD= 0.27497 RMSD= 0.43539
    case wb97:
      s8 = 7.11022468;
      a1 = 0.76423345;
      a2 = 8.44559334;
      break;
      // Fitset: MD= -0.12818 MAD= 0.36163 RMSD= 0.50023
    case wb97x:
      s8 = 0.38815338;
      a1 = 0.47448629;
      a2 = 6.91367384;
      break;
      // Fitset: MD= -0.20303 MAD= 0.34880 RMSD= 0.54190
    case x3lyp:
      s8 = 1.55067492;
      a1 = 0.19818545;
      a2 = 5.61262748;
      break;
      // Fitset: MD= -0.17267 MAD= 0.32089 RMSD= 0.51125
    case xlyp:
      s8 = 1.51577878;
      a1 = 0.10026585;
      a2 = 5.37506460;
      break;
      // Fitset: MD= -0.06152 MAD= 0.27318 RMSD= 0.38360
  }

  par.s6 = s6;
  par.s8 = s8;
  par.s9 = s9;
  par.s10 = s10;
  par.a1 = a1;
  par.a2 = a2;
  par.alp = alp;

  return par;
}

int d4par(
  const std::string func,
  dftd::dparam &par,
  const bool latm/* = true*/
) {
  auto num = get_dfunc(func);
  if (latm) {
    par = get_d4eeqbjatm_2019_parameter(num);
  } else {
    par = get_d4eeqbjmbd_2019_parameter(num);
  }
  return EXIT_SUCCESS;
}

}  // namespace dftd
