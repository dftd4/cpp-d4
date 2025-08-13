
#ifndef TEST_MULTICHARGE_H
#define TEST_MULTICHARGE_H

static const double param_eeqbc_chi_ref[2]{
  1.7500687479,
  0.6679573735
};

static const double param_eeqbc_rvdw_ref[3]{
  2.1823,
  4.3286,
  4.1429
}; 

static const double qvec_reference[]{
     0.4757830909, -0.0426540501, -0.3778712260, -0.0967376090,
    -0.1733641170,  0.1086601010, -0.1136284484, -0.3179396996,
    -0.2456555247,  0.1761065724,  0.1145108507, -0.1222410255,
    -0.0144595425,  0.2577820828, -0.1117775795,  0.4834861246
  };

extern int test_multi();

#endif  // TEST_MULTICHARGE_H