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

tests = [
  'disp2',
  'disp',
  'ghost',
  'grad',
  'ncoord',
  'param',
]

test_srcs = files(
  'main.cpp',
  'molecules.cpp',
  'util.cpp',
)
foreach t : tests
  test_srcs += files('test_@0@.cpp'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [cpp_d4_dep],
  link_language: 'cpp',
)

foreach t : tests
  test(t, tester, args: t)
endforeach
