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
#ifndef MOLECULES_H
#define MOLECULES_H

// MB16_43: 01
static const int mb16_43_01_charge{0};
static const char mb16_43_01_atoms[16][3] {
  "Na", "H", "O", "H", "F", "H", "H", "O", "N", "H", "H", "Cl", "B",
  "B", "N", "Al",
};
static const double mb16_43_01_coord[16*3] {
  -1.85528263484662, +3.58670515364616, -2.41763729306344,
  +4.40178023537845, +0.02338844412653, -4.95457749372945,
  -2.98706033463438, +4.76252065456814, +1.27043301573532,
  +0.79980886075526, +1.41103455609189, -5.04655321620119,
  -4.20647469409936, +1.84275767548460, +4.55038084858449,
  -3.54356121843970, -3.18835665176557, +1.46240021785588,
  +2.70032160109941, +1.06818452504054, -1.73234650374438,
  +3.73114088824361, -2.07001543363453, +2.23160937604731,
  -1.75306819230397, +0.35951417150421, +1.05323406177129,
  +5.41755788583825, -1.57881830078929, +1.75394002750038,
  -2.23462868255966, -2.13856505054269, +4.10922285746451,
  +1.01565866207568, -3.21952154552768, -3.36050963020778,
  +2.42119255723593, +0.26626435093114, -3.91862474360560,
  -3.02526098819107, +2.53667889095925, +2.31664984740423,
  -2.00438948664892, -2.29235136977220, +2.19782807357059,
  +1.12226554109716, -1.36942007032045, +0.48455055461782,
};

// ROST61: m1
static const int rost61_m1_charge{0};
static const char rost61_m1_atoms[22][3] {
  "c", "c", "c", "h", "c", "c", "h", "h", "h", "h", "c", "c", "c", "h",
  "c", "c", "h", "h", "h", "h", "ti", "h",
}; 
static const double rost61_m1_coord[22*3] {
  -1.71626628838550,  1.81517342801759,  1.07097099223840,
  -2.51566406323399,  3.03001238435168, -1.17663836146950,
  -3.60349996064212,  1.18149002455176, -2.75730692172184,
  -2.32680525082228,  5.01060073821355, -1.60490109782415,
  -3.46303952033678, -1.16020920405222, -1.49616837438297,
  -2.29891779482473, -0.78713713803762,  0.87106590797531,
  -1.92987748807901, -2.21384610996610,  2.27511460260479,
  -0.83020368961167,  2.71965222693500,  2.66216051001990,
  -4.35678220830401,  1.48747706230290, -4.61847148201756,
  -4.04788499157012, -2.95326796832416, -2.26881632437650,
   4.64052109664389,  1.78126853688051, -3.77294172207966,
   4.18011026520508,  2.98511706707145, -1.45305517088802,
   4.13214560211617,  1.11344011809522,  0.45760731483005,
   3.85989860827593,  4.97729166445008, -1.18037453164844,
   4.56810500061798, -1.26128365002214, -0.68598818680096,
   4.86832535633142, -0.83471377648079, -3.30933738716750,
   5.13998542776287, -2.26171560296739, -4.73586884461798,
   4.70588989810968,  2.67256562769832, -5.59852359075679,
   3.82913510402041,  1.43996458655102,  2.44202982151600,
   4.68804786915785, -3.05243319181537,  0.27393641783910,
   0.73073048105294,  0.23022005138048, -2.18757235919326,
   0.23556166767074,  0.10105333453460, -5.42771553810106,
};

#endif // MOLECULES_H
