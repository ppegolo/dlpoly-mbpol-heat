#ifndef POLY_2B_A1_A1_V1X_H
#define POLY_2B_A1_A1_V1X_H

namespace A1_A1 {

//
// this is the polynomial used by x2b_h2o_ion_v1<5> (including gradients)
//

struct poly_2b_A1_A1_v1x {
    static const unsigned degree = 13;
    static const unsigned n_vars = 1;
    static const unsigned size = 13;

    static double eval(const double a[13],
                       const double x[1]);

    static double eval(const double a[13],
                       const double x[1],
                             double g[1]);
};

} // namespace h2o_ion

#endif // POLY_2B_H2O_ION_V1X_H
