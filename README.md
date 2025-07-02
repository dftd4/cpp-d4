# C++ Port of DFT-D4

[![License](https://img.shields.io/github/license/dftd4/cpp-d4)](https://github.com/dftd4/cpp-d4/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/dftd4/cpp-d4)](https://github.com/dftd4/cpp-d4/releases/latest)

This project is a port of the [`dftd4`](https://github.com/dftd4/dftd4) project
to C++ and provides the D4(EEQ)-ATM method.

## Building This Project

This project is build with `meson`, to setup and perform a build run:

```bash
meson setup _build
meson compile -C _build
```

Run the test suite with:

```bash
meson test -C _build --print-errorlogs
```

## Citations

- E. Caldeweyher, C. Bannwarth and S. Grimme, _J. Chem. Phys._, **2017**, 147, 034112. DOI: [10.1063/1.4993215](https://dx.doi.org/10.1063/1.4993215)
- E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth and S. Grimme, _J. Chem. Phys._, **2019**, 150, 154122. DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
- E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, _Phys. Chem. Chem. Phys._, **2020**, 22, 8499-8512. DOI: [10.1039/D0CP00502A](https://doi.org/10.1039/D0CP00502A)

## License

### DFT-D4 Code

DFT-D4 is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DFT-D4 is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
GNU Lesser General Public License for more details.

### ORCA-derived Code in `dftd_cblas.h` and `dftd_matrix.h`

The files [`include/dftd_cblas.h`](https://github.com/dftd4/cpp-d4/blob/main/include/dftd_cblas.h) and
[`include/dftd_matrix.h`](https://github.com/dftd4/cpp-d4/blob/main/include/dftd_matrix.h)
contain code adapted from the ORCA quantum chemistry program,
which is developed by the group of Prof. Frank Neese at the Max-Planck-Institut für Kohlenforschung,
Mülheim an der Ruhr and FAccTs GmbH. ORCA is licensed by the
Max-Planck-Institut für Kohlenforschung and FAccTs GmbH.

The inclusion of ORCA code in these files is done with the explicit permission
of the ORCA developers. This code remains subject to the licensing terms
of ORCA, which allow free academic use but require a separate license for
industrial or commercial use.

For reuse or licensing of the code in `dftd_cblas.h` and `dftd_matrix.h`, please contact
the ORCA team at the Max-Planck-Institut für Kohlenforschung (https://orcaforum.kofo.mpg.de/)
or FAccTs GmbH. (https://www.faccts.de/).
