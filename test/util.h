#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <dftd_geometry.h>


extern int get_molecule(int n, const char atoms[][3], const double coord[], dftd::TMolecule& mol);

extern bool check(double actual, double expected, double epsilon = 1e-12, bool rel = false);
extern bool check(float actual, float expected, float epsilon = 1e-6, bool rel = false);

extern int element(const std::string& sym);

#endif // UTIL_H