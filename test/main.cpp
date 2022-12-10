#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>

#include "test_disp.h"
#include "test_ncoord.h"


enum test {
  invalid,
  disp,
  ncoord,
};

test get_tests(std::string name) {
  static const std::map<std::string, test> testStrings {
      {"disp", disp},
      {"ncoord", ncoord},
  };

  std::string test = name;
  transform(name.begin(), name.end(), test.begin(), ::tolower);
  auto iter = testStrings.find(test);
  if (iter != testStrings.end()) {
    return iter->second;
  }
  return invalid;
};


int main(int argc, char *argv[]) {
  switch (get_tests(argv[1])) {
    default:
      return EXIT_FAILURE;
    case ncoord:
      return test_ncoord(); 
    case disp:
      return test_disp(); 
  }
}