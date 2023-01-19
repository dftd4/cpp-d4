# C++ Port of DFT-D4

[![License](https://img.shields.io/github/license/dftd4/cpp-d4)](https://github.com/dftd4/cpp-d4/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/dftd4/cpp-d4)](https://github.com/dftd4/cpp-d4/releases/latest)

This project is a port of the [`dftd4`](https://github.com/dftd4/dftd4) project
to C++ and provides the D4(EEQ)-ATM method.

## Requirements
cpp-d4 depends on cblas and lapacke, the C implementation of blas and lapack.


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

In addition to `meson`, one can use CMake to build cpp-d4 and include it into external projects (see for example [https://github.com/conradhuebler/curcuma](https://github.com/conradhuebler/curcuma)). 

```bash
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make
```
The tests can be started with

```bash
make test
```
In case, cblas and lapacke can not be found using cmake (which for example happens on Ubuntu), one can specify one include and library path for both dependencies.

```bash
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. -DD4LIBS_DIR=/usr/lib/ -DD4INCLUDE_DIR=/usr/include/
```

```bash
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. -DD4LIBS_DIR=$(find /usr -name 'cblas.so' |head -n1 |sed 's/cblas.so//g') -DD4INCLUDE_DIR=$(find /usr -name 'cblas.h' |head -n1 |sed 's/cblas.h//g')
```

Using CMake it is then easy to include cpp-d4 in an own project.

```cmake
if(NOT DEFINED D4LIBS_DIR AND NOT DEFINED D4INCLUDE_DIR)
  	find_package(LAPACKE REQUIRED)
	find_package(CBLAS CONFIG REQUIRED)
else()
	include_directories(${INCLUDE_DIR})
endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)

add_subdirectory(cpp-d4)

if(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
    target_link_libraries(yourproject PUBLIC libcpp_d4 cblas lapacke )
else()
    target_link_libraries(yourproject PUBLIC libcpp_d4 ${LAPACKE_LIBRARIES}  ${CBLAS_LIBRARIES} )
endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
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
