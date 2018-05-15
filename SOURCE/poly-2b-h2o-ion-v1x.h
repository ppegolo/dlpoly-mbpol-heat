#ifndef POLY_2B_H2O_ION_V1X_H
#define POLY_2B_H2O_ION_V1X_H

namespace h2o_ion {

//
// this is the polynomial used by x2b_h2o_ion_v1<5> (including gradients)
//

struct poly_2b_h2o_ion_v1x {
    static const unsigned degree = 5;
    static const unsigned n_vars = 8;
    static const unsigned size = 429;

    static double eval(const double a[428],
                       const double x[8]);

    static double eval(const double a[429],
                       const double x[8],
                             double g[8]);
};

} // namespace h2o_ion

#endif // POLY_2B_H2O_ION_V1X_H
