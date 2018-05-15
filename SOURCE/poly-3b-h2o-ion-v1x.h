#ifndef POLY_3B_H2O_ION_V1X_H
#define POLY_3B_H2O_ION_V1X_H

namespace h2o_ion {

//
// this is the polynomial used by x3b_h2o_ion_v2 (including gradients)
//

struct poly_3b_h2o_ion_v1x {
    static const unsigned n_vars = 21;
    static const unsigned size = 1831;

    static double eval(const double a[1831],
                       const double x[21]);

    static double eval(const double a[1831],
                       const double x[21],
                             double g[21]);
};

} // namespace h2o_ion

#endif // POLY_3B_H2O_ION_V2X_H
