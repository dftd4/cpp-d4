# C++ Port of DFT-D4 (for ORCA)

[![License](https://img.shields.io/github/license/dftd4/cpp-d4)](https://github.com/dftd4/cpp-d4/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/dftd4/cpp-d4)](https://github.com/dftd4/cpp-d4/releases/latest)

This project is a port of the [`dftd4`](https://github.com/dftd4/dftd4) project
to C++ and provides the D4(EEQ)-ATM method.

**NOTE:** This branch contains some adaptations for compatibility with the ORCA Quantum Chemistry program. This inludes some renaming of variables and functions, and, most notably, the handling of ghost atoms.

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

DFT-D4 is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DFT-D4 is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
GNU Lesser General Public License for more details.
