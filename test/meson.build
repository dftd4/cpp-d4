tests = [
  'disp',
  'ncoord',
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
