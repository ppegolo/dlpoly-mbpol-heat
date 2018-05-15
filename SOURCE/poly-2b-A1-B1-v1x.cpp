
#include "poly-2b-A1-B1-v1x.h"

//namespace mb_system_fit {

namespace A1_B1 {

double poly_2b_A1_B1_v1x::eval(const double a[13], const double x[1],
                               double g[1])
{
    const double t1 = a[0];
    const double t2 = a[1];
    const double t3 = a[2];
    const double t4 = a[3];
    const double t5 = a[4];
    const double t6 = a[5];
    const double t7 = a[6];
    const double t8 = a[7];
    const double t9 = a[8];
    const double t10 = a[9];
    const double t11 = a[10];
    const double t15 = x[0];
    const double t13 = a[12]*t15;
    const double t14 = a[11];
    const double t16 = (t13+t14)*t15;
    const double t18 = (t11+t16)*t15;
    const double t20 = (t10+t18)*t15;
    const double t22 = (t9+t20)*t15;
    const double t24 = (t8+t22)*t15;
    const double t26 = (t7+t24)*t15;
    const double t28 = (t6+t26)*t15;
    const double t30 = (t5+t28)*t15;
    const double t32 = (t4+t30)*t15;
    const double t34 = (t3+t32)*t15;
    const double t36 = (t2+t34)*t15;
    g[0] = (((((((((((2.0*t13+t14)*t15+t11+t16)*t15+t10+t18)*t15+t9+t20)
*t15+t8+t22)*t15+t7+t24)*t15+t6+t26)*t15+t5+t28)*t15+t4+t30)*t15+t3+t32)*t15+t2
+t34)*t15+t1+t36;
    return (t1+t36)*t15;
}
} // namespace mb_system
