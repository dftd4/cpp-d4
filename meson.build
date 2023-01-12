# This file is part of cpp-d4.
#
# Copyright (C) 2019 Sebastian Ehlert, Marvin Friede
#
# cpp-d4 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cpp-d4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cpp-d4.  If not, see <https://www.gnu.org/licenses/>.

project(
  'cpp-d4',
  'cpp',
  version: '2.0',
  license: 'LGPL3-or-later',
  meson_version: '>=0.45',
)

cxx = meson.get_compiler('cpp')

deps = []
foreach lib: ['cblas', 'lapacke']
  this_dep = dependency(lib, required: false)
  if not this_dep.found()
    this_dep = cxx.find_library(lib)
  endif
  deps += this_dep
endforeach

srcs = files(
  'src/dftd_cutoff.cpp',
  'src/dftd_damping.cpp',
  'src/dftd_dispersion.cpp',
  'src/dftd_eeq.cpp',
  'src/dftd_model.cpp',
  'src/dftd_ncoord.cpp',
  'src/damping/atm.cpp',
  'src/damping/rational.cpp',
)

cpp_d4_inc = include_directories('include')

cpp_d4_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  include_directories: cpp_d4_inc,
  dependencies: deps,
  install: true,
)

install_headers(
  'include/dftd_cutoff.h',
  'include/dftd_damping.h',
  'include/dftd_dispersion.h',
  'include/dftd_eeq.h',
  'include/dftd_geometry.h',
  'include/dftd_matrix.h',
  'include/dftd_model.h',
  'include/dftd_ncoord.h',
  'include/dftd_parameters.h',
  'include/damping/atm.h',
  'include/damping/rational.h',
)

cpp_d4_dep = declare_dependency(
  link_with: cpp_d4_lib,
  include_directories: cpp_d4_inc,
  dependencies: deps,
)

executable(
  meson.project_name(),
  sources: files(
    'src/program_dftd.cpp',
    'src/dftd_readxyz.cpp',
  ),
  dependencies: cpp_d4_dep,
  install: true,
)

subdir('test')
