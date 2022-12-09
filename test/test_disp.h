#ifndef TEST_DISP_H
#define TEST_DISP_H

static const double mb16_43_01_ref_energy {
  -2.5882588037023E-02
};

static const double rost61_m1_ref_energy {
  -3.4287391104745E-02
};

extern int test_energy(
  const int n,
  const char atoms[][3],
  const double coord[],
  const double ref_cn
);

extern int test_disp(void);

#endif // TEST_DISP_H
