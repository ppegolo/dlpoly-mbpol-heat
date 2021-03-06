#include <cmath>
#include <algorithm>

#include "poly-2b-h2o-ion-v1x.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

struct variable {
    double v_exp(const double& r0, const double& k,
                 const double* xcrd, int o, int x );

    double v_coul(const double& r0, const double& k,
                  const double* xcrd, int o, int x );

    void grads(const double& gg, double* xgrd, int o, int x ) const;

    double g[3]; // diff(value, p1 - p2)
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_exp(const double& r0, const double& k,
                       const double* xcrd, int o, int x )
{
    g[0] = xcrd[o++] - xcrd[x++];
    g[1] = xcrd[o++] - xcrd[x++];
    g[2] = xcrd[o]   - xcrd[x];

    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);

    const double exp1 = std::exp(k*(r0 - r));
    const double gg = - k*exp1/r;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return exp1;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_coul(const double& r0, const double& k,
                        const double* xcrd, int o, int x)
{
    g[0] = xcrd[o++] - xcrd[x++];
    g[1] = xcrd[o++] - xcrd[x++];
    g[2] = xcrd[o]   - xcrd[x];

    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    const double r = std::sqrt(rsq);

    const double exp1 = std::exp(k*(r0 - r));
    const double rinv = 1.0/r;
    const double val = exp1*rinv;

    const double gg = - (k + rinv)*val*rinv;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return val;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void variable::grads(const double& gg, double* xgrd, int o, int x) const
{
    for (int i = 0; i < 3; ++i) {
        const double d = gg*g[i];

        xgrd[o++] += d;
        xgrd[x++] -= d;
    }
}

//----------------------------------------------------------------------------//

struct monomer {
    double oh1[3];
    double oh2[3];

    void setup(const double* ohh,
               const double& in_plane_g, const double& out_of_plane_g,
               double x1[3], double x2[3]);

    void grads(const double* g1, const double* g2,
               const double& in_plane_g, const double& out_of_plane_g,
               double* grd) const;
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::setup(const double* ohh,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* x1, double* x2)
{
    for (int i = 0; i < 3; ++i) {
        oh1[i] = ohh[i + 3] - ohh[i];
        oh2[i] = ohh[i + 6] - ohh[i];
    }

    const double v[3] = {
        oh1[1]*oh2[2] - oh1[2]*oh2[1],
        oh1[2]*oh2[0] - oh1[0]*oh2[2],
        oh1[0]*oh2[1] - oh1[1]*oh2[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double in_plane = ohh[i] + 0.5*in_plane_g*(oh1[i] + oh2[i]);
        const double out_of_plane = out_of_plane_g*v[i];

        x1[i] = in_plane + out_of_plane;
        x2[i] = in_plane - out_of_plane;
    }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::grads(const double* g1, const double* g2,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* grd) const
{
    const double gm[3] = {
        g1[0] - g2[0],
        g1[1] - g2[1],
        g1[2] - g2[2]
    };

    const double t1[3] = {
        oh2[1]*gm[2] - oh2[2]*gm[1],
        oh2[2]*gm[0] - oh2[0]*gm[2],
        oh2[0]*gm[1] - oh2[1]*gm[0]
    };

    const double t2[3] = {
        oh1[1]*gm[2] - oh1[2]*gm[1],
        oh1[2]*gm[0] - oh1[0]*gm[2],
        oh1[0]*gm[1] - oh1[1]*gm[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double gsum = g1[i] + g2[i];
        const double in_plane = 0.5*in_plane_g*gsum;

        const double gh1 = in_plane + out_of_plane_g*t1[i];
        const double gh2 = in_plane - out_of_plane_g*t2[i];

        grd[i + 0] += gsum - (gh1 + gh2); // O
        grd[i + 3] += gh1; // H1
        grd[i + 6] += gh2; // H2
    }
}

////////////////////////////////////////////////////////////////////////////////

struct x2b_h2o_ion_v1x_p {
    x2b_h2o_ion_v1x_p();

    typedef h2o_ion::poly_2b_h2o_ion_v1x poly_type;

    double operator()(const double* w, const double* x,
                      double* g1, double* g2) const;

protected:
    double m_k_HH_intra;
    double m_k_OH_intra;

    double m_k_XH_coul;
    double m_k_XO_coul;

    double m_k_XLp_main;

    double m_d_HH_intra;
    double m_d_OH_intra;

    double m_d_XH_coul;
    double m_d_XO_coul;

    double m_d_XLp_main;

    double m_in_plane_gamma;
    double m_out_of_plane_gamma;

public:
    double m_r2i;
    double m_r2f;

    double f_switch(const double&, double&) const; // X-O separation

protected:
    double m_poly[poly_type::size];
};

////////////////////////////////////////////////////////////////////////////////

double x2b_h2o_ion_v1x_p::f_switch(const double& r, double& g) const
{
    if (r > m_r2f) {
        g = 0.0;
        return 0.0;
    } else if (r > m_r2i) {
        const double t1 = M_PI/(m_r2f - m_r2i);
        const double x = (r - m_r2i)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

//----------------------------------------------------------------------------//

double x2b_h2o_ion_v1x_p::operator()
    (const double* w, const double* x, double* g1, double* g2) const
{
    // the switch

    const double dXO[3] = {w[0] - x[0],
                           w[1] - x[1],
                           w[2] - x[2]};

    const double rXOsq = dXO[0]*dXO[0] + dXO[1]*dXO[1] + dXO[2]*dXO[2];
    const double rXO = std::sqrt(rXOsq);

    if (rXO > m_r2f)
        return 0.0;

    // offsets

    const int O  = 0;
    const int H1 = 3;
    const int H2 = 6;

    const int X   = 9;
    
    const int Lp1 = 12;
    const int Lp2 = 15;

    double xcrd[18]; // coordinates including extra-points

    std::copy(w, w + 9, xcrd);
    std::copy(x , x  + 9, xcrd + 9);

    // the extra-points

    monomer ma;

    ma.setup(xcrd + O,
             m_in_plane_gamma, m_out_of_plane_gamma,
             xcrd + Lp1, xcrd + Lp2);

    // variables

//    const double d0_intra = 1.0;  //TODO MBpol values 
//    const double d0_inter = 4.0;

    double v[ 8]; // stored separately (gets passed to poly::eval)

    variable ctxt[8 ];

    v[0] = ctxt[0].v_exp(m_d_HH_intra, m_k_HH_intra, xcrd, H1, H2);

    v[1] = ctxt[1].v_exp(m_d_OH_intra, m_k_OH_intra, xcrd, O, H1);
    v[2] = ctxt[2].v_exp(m_d_OH_intra, m_k_OH_intra, xcrd, O, H2);

    v[3] = ctxt[3].v_coul(m_d_XH_coul, m_k_XH_coul, xcrd, X, H1);
    v[4] = ctxt[4].v_coul(m_d_XH_coul, m_k_XH_coul, xcrd, X, H2);
                            
    v[5] = ctxt[5].v_coul(m_d_XO_coul, m_k_XO_coul, xcrd, X, O);

    v[6] = ctxt[6].v_exp(m_d_XLp_main, m_k_XLp_main, xcrd, X, Lp1);
    v[7] = ctxt[7].v_exp(m_d_XLp_main, m_k_XLp_main, xcrd, X, Lp2);

    double g[8 ];
    const double E_poly = h2o_ion::poly_2b_h2o_ion_v1x::eval(m_poly, v, g);

    double xgrd[18];
    std::fill(xgrd, xgrd + 18, 0.0);

   ctxt[0].grads(g[0], xgrd, H1, H2);
                       
    ctxt[1].grads(g[1], xgrd, O, H1);
    ctxt[2].grads(g[2], xgrd, O, H2);
                       
    ctxt[3].grads(g[3], xgrd, X, H1);
    ctxt[4].grads(g[4], xgrd, X, H2);

    ctxt[5].grads(g[5], xgrd, X, O);

    ctxt[6].grads(g[6], xgrd, X, Lp1);
    ctxt[7].grads(g[7], xgrd, X, Lp2);

    // distribute gradients w.r.t. the X-points

    ma.grads(xgrd + Lp1, xgrd + Lp2,
             m_in_plane_gamma, m_out_of_plane_gamma,
             xgrd + O);

    // the switch

    double gsw;
    const double sw = f_switch(rXO, gsw);

    for (int i = 0; i < 9; ++i) {
        g1[i] = sw*xgrd[i];
    }

    for (int i = 0; i < 3; ++i) {
        g2[i] += sw*xgrd[i + 9];
    }

    // gradient of the switch

    gsw *= E_poly/rXO;
    for (int i = 0; i < 3; ++i) {
        const double d = gsw*dXO[i];
        g1[i] += d;
        g2[i] -= d;
    }

    return sw*E_poly;
}

//----------------------------------------------------------------------------//

//fit-fullpolargrid-fixedwaterparams  12/20/16
static const double the_poly[] = {
 3.172252449777494e+02, // 0
 1.534926181980350e+02, // 1
-2.507456307693576e+01, // 2
-1.300880086477685e+02, // 3
 4.077107977837585e+02, // 4
-2.426747550116293e+01, // 5
 1.178105617614842e+02, // 6
-2.509574395977389e+02, // 7
 4.444360706737240e+01, // 8
-1.752548109046285e+01, // 9
 1.155678093539773e+01, // 10
-1.535599032358893e+02, // 11
-6.162309916155965e+02, // 12
-1.551519174619939e+03, // 13
 5.316443296119055e+02, // 14
 3.811406989049416e+01, // 15
 8.517784295925232e+01, // 16
 8.150766110975600e+01, // 17
 3.958344649485545e+02, // 18
 2.111267854685706e+02, // 19
-2.675324384876650e+01, // 20
 2.912063456929841e+02, // 21
-1.257158281852784e+02, // 22
 8.196210695093701e+01, // 23
-2.296326739422743e+00, // 24
 3.392475749725058e+02, // 25
-2.166766051326271e+01, // 26
 5.428107048397512e+00, // 27
-4.312323469153977e+00, // 28
-4.895305818116886e+01, // 29
-2.228793016216545e+02, // 30
-1.146755491725045e+01, // 31
 5.268349483940045e+01, // 32
 8.502233791583340e+01, // 33
 8.315699969804744e+01, // 34
-4.139778494167293e+00, // 35
 1.948491231381103e+02, // 36
-3.646236411053430e+01, // 37
 5.796169234802060e+00, // 38
 8.756697162881484e+01, // 39
-3.737857414104710e+01, // 40
-1.115475640776338e+03, // 41
-2.403319750263057e-01, // 42
 1.429852265237522e+02, // 43
-2.693341704814821e+02, // 44
-6.283449608193100e+00, // 45
 3.242098079295404e+02, // 46
-1.075835190365443e+02, // 47
 2.822040693330885e+02, // 48
 3.798411499309777e+02, // 49
 5.073024288785667e+01, // 50
 7.214692962627250e+01, // 51
-2.891631306442829e+02, // 52
-2.200937173211147e+02, // 53
-4.781381442551095e+02, // 54
 3.343753107645887e+02, // 55
 1.322088028081601e+00, // 56
 9.460230985694152e+01, // 57
 3.258374315616745e+02, // 58
 7.877815372177309e+01, // 59
 6.816193554638027e+02, // 60
 2.016423328822777e+02, // 61
 4.792713104023705e+02, // 62
-3.218679521163862e+02, // 63
 8.310732074665665e+00, // 64
-3.106252068250761e+02, // 65
-7.710082041957278e+01, // 66
-1.903760711865668e+01, // 67
 2.112269041748477e+00, // 68
-7.246140066432903e+01, // 69
-2.024251165409403e-01, // 70
 3.175926161964438e+01, // 71
-5.919644047752028e+00, // 72
 1.529440787482422e+02, // 73
 1.107767562539906e+02, // 74
 1.319252637161772e+01, // 75
-3.896036827676173e-02, // 76
 3.192168824772019e+02, // 77
-5.397646281079684e+02, // 78
 1.040632548038264e+03, // 79
 2.879005218444141e+01, // 80
 1.358274270255133e+01, // 81
 2.047598113789306e+02, // 82
-6.064507705245326e-01, // 83
-8.313252402740672e+01, // 84
-5.937293442828944e+00, // 85
-2.299979954237862e+01, // 86
-3.598883335852660e+02, // 87
-5.367071543458617e+01, // 88
-4.748000437324333e+01, // 89
-5.665742982010274e+01, // 90
 3.679688651744412e+00, // 91
-6.562886550052022e-01, // 92
-2.114683815977760e+02, // 93
 6.063090320032241e+01, // 94
-9.227818769268108e+02, // 95
-2.420964827162172e+02, // 96
 2.650641553847398e+01, // 97
-1.955482163682205e+01, // 98
-8.735821081925626e+01, // 99
 5.040304490004117e+02, // 100
 7.189459130933605e+01, // 101
-3.188196666712204e+01, // 102
-6.405258526203494e-01, // 103
-4.130823484847314e+02, // 104
 3.102776191006105e+01, // 105
 5.217207629596037e+01, // 106
 2.351559956973298e+01, // 107
-2.093556269863195e+02, // 108
 2.192712020749362e+02, // 109
-7.943515182573836e+00, // 110
 9.985407732886348e+01, // 111
 4.220185177538920e+00, // 112
 3.549226096753247e+01, // 113
 1.156027364597560e+02, // 114
-1.419445625339465e+00, // 115
 1.690412768645464e+01, // 116
-2.035827486547527e+01, // 117
 3.809257247651406e+02, // 118
 1.859689852604427e+01, // 119
-8.289287293396113e+00, // 120
 1.018205885642080e+01, // 121
-4.315212317892971e-02, // 122
 5.946020253118377e+00, // 123
 5.345333959835094e+00, // 124
-2.406740053250880e+00, // 125
-1.153177982672368e+00, // 126
 2.417966883697390e+02, // 127
-1.194463250361711e+00, // 128
-1.366667171329481e+03, // 129
-1.744573620806214e+00, // 130
-6.411570033034766e+01, // 131
 1.285226467886173e+00, // 132
-3.910961559639568e+02, // 133
-1.440988694665527e+01, // 134
-1.268403640189199e+01, // 135
 1.783165264135599e+02, // 136
 5.531146124932191e+02, // 137
 4.402881897397054e+01, // 138
-5.843631032537530e+01, // 139
-1.530564834236027e+01, // 140
 4.540146888328814e+00, // 141
 5.619433647429610e+00, // 142
 2.042101329213225e+02, // 143
 4.769377301080399e+01, // 144
 9.814084140708697e-01, // 145
-1.019290522989328e+01, // 146
-6.171875699866097e+01, // 147
-1.552963713021786e+02, // 148
-1.834646777761779e-01, // 149
-2.073143884040286e+01, // 150
-2.212850594565412e+00, // 151
 1.534238190178292e+01, // 152
-2.356298619913602e+00, // 153
 1.115538403649581e-02, // 154
-2.346622715299085e+00, // 155
-9.958955597173723e+01, // 156
 8.401989673467900e+02, // 157
-6.465607982237827e-03, // 158
-5.555284785555767e+01, // 159
 1.723311187545044e+00, // 160
 2.200292169237219e+02, // 161
-6.083841664237525e+02, // 162
 9.909144591381281e+02, // 163
 3.287473687001403e+01, // 164
 4.005897177509266e+01, // 165
 3.172254625781076e+02, // 166
-1.469959561184984e+00, // 167
 9.226326520642263e+01, // 168
 1.056181185826668e+03, // 169
-8.201802790795174e+00, // 170
-2.643226498941732e+02, // 171
-2.616082004802192e+02, // 172
-5.721549060086209e+01, // 173
 5.535069295961729e+00, // 174
-3.452316283873698e+01, // 175
 1.129697160677676e+01, // 176
-4.686380667760980e+01, // 177
-6.091721797413991e+01, // 178
-4.187017190779328e+02, // 179
 1.145803711070823e-01, // 180
 4.188263717229673e+01, // 181
-7.123651332523999e+00, // 182
-8.565914976842839e+01, // 183
-4.702862526424481e+00, // 184
 6.176680888294161e+02, // 185
 3.373395290222005e+00, // 186
-1.176786595868546e+01, // 187
 6.808512056506189e+00, // 188
 5.126527758111734e-02, // 189
 1.756833629342558e+01, // 190
 1.642943568389019e+01, // 191
-7.604226037262050e-02, // 192
-5.584641478328830e-01, // 193
-1.515958873340970e+02, // 194
-2.718716903580688e+02, // 195
-1.465222236858225e-02, // 196
-2.363136484568214e+02, // 197
-5.803736961770692e+01, // 198
 1.761091589925620e+00, // 199
 6.188317771359003e+02, // 200
-3.461531646690386e-01, // 201
-3.051637676170656e+00, // 202
-3.625113391591061e-01, // 203
 3.918109049224480e+00, // 204
 2.084410480732735e+01, // 205
 5.620878561951391e+00, // 206
-7.726779682910683e-01, // 207
-2.650829549183901e+00, // 208
-3.646395961600407e+01, // 209
-2.750356876073431e-04, // 210
 2.109157315652825e-03, // 211
 2.703744502719805e+02, // 212
 1.120115547504763e-02, // 213
 1.255606997384222e+00, // 214
 1.308814649164288e+00, // 215
-8.341705968404458e+00, // 216
 2.158902246134387e+00, // 217
 6.400583764160610e+00, // 218
-1.944239366962265e+00, // 219
 7.529302585347303e+00, // 220
-6.253274577457094e-02, // 221
 1.583644515067475e-01, // 222
 8.454206633419482e-02, // 223
-2.925830204746839e-01, // 224
 5.918788946963590e+02, // 225
-9.595059427673077e-01, // 226
-3.255234188473422e-01, // 227
 6.931712711718477e+00, // 228
 2.150405370140902e+01, // 229
-6.253856373853507e+00, // 230
 5.186024489190666e-02, // 231
 2.173427341004596e+01, // 232
-3.114045626795448e+00, // 233
-1.679354413330871e+00, // 234
 7.183509831302876e+00, // 235
-1.361622560854541e-03, // 236
-1.048070436825830e-01, // 237
-1.479764474695675e+01, // 238
-5.320665103431333e+00, // 239
 4.644039331012360e+01, // 240
-8.549700573985191e-01, // 241
-1.589873520081237e+00, // 242
-3.044201751454393e+01, // 243
 2.389297550618046e+00, // 244
-2.084922572604866e-02, // 245
-4.170140249851713e-01, // 246
-4.660842459398098e-01, // 247
 5.925633823692532e+01, // 248
 1.014259492581146e+00, // 249
-3.329082834385294e+00, // 250
 8.406442332037650e-02, // 251
-6.263867301994875e+01, // 252
 4.510497075860565e+01, // 253
 3.804527113625460e+00, // 254
 3.875723264723032e-01, // 255
-1.095358934391015e+02, // 256
-2.196601298464055e+00, // 257
 2.308352373995986e-02, // 258
-3.295768920482142e+00, // 259
-1.829598133276994e+01, // 260
 1.069033784256045e-02, // 261
-1.536838911633021e-01, // 262
 8.581369676367044e-01, // 263
-1.902329694258683e-01, // 264
 3.845489283309411e+01, // 265
-6.764852031610069e+01, // 266
 3.779618229799145e+01, // 267
 6.603126136060914e-04, // 268
-3.609542281052842e-05, // 269
-1.388801573744304e+01, // 270
 5.228337678863431e-03, // 271
-1.111782947965339e+00, // 272
-2.348844099827950e-01, // 273
 2.682344877235188e-02, // 274
 3.784192733517906e-02, // 275
-8.492616865102207e+02, // 276
 3.434651352275791e-05, // 277
 6.060237791029212e+00, // 278
 7.607073768299605e+02, // 279
-9.455163647375845e+01, // 280
-1.140713152596373e+01, // 281
-1.499882859186185e+01, // 282
-8.133533754799494e+00, // 283
 5.904016804695048e+00, // 284
-1.819507799956342e-02, // 285
 2.532620988353231e+01, // 286
 2.138728530012947e-01, // 287
-3.593795940263081e-03, // 288
 1.126785371924476e+02, // 289
-6.688979268745762e-01, // 290
 3.451651325841812e-01, // 291
-1.407771603640340e-01, // 292
 1.050012874000041e-01, // 293
 2.363179879070822e+01, // 294
-2.474735716390669e+02, // 295
 2.965674061399131e+00, // 296
 2.975552964356542e+01, // 297
 5.154758590337605e-02, // 298
-2.530254756490458e-04, // 299
 3.324008160553047e-01, // 300
-2.167189339524576e+02, // 301
-2.079384725124306e+02, // 302
-8.612651774181651e-01, // 303
 3.221792280249202e-01, // 304
-2.773193517599982e+02, // 305
-4.399551099490152e+01, // 306
-3.689370125560649e-01, // 307
 2.791736358342926e+02, // 308
 6.413636786354430e+00, // 309
-6.380786403670422e-01, // 310
 5.626357478574155e-02, // 311
-1.350438774470658e-02, // 312
 4.570509221648325e-01, // 313
-7.926175908845177e-03, // 314
-2.872709746812037e+00, // 315
 1.310638329366798e+01, // 316
-5.010624507153254e-01, // 317
-1.488043910081159e+01, // 318
 2.576474795922340e+01, // 319
 8.260014060876823e-02, // 320
 4.896242417230783e-01, // 321
-3.835260875620855e+01, // 322
 5.317385557703619e-01, // 323
 3.550517901963035e-01, // 324
 6.417079835986557e+00, // 325
-8.692713015149829e+00, // 326
 3.986090535410609e+00, // 327
 4.033669915239619e+00, // 328
 1.571510652712739e+00, // 329
-3.920278083534410e-04, // 330
-1.570474501262159e+01, // 331
 2.341860187396334e+01, // 332
 5.563143046521513e-01, // 333
-4.127551266716997e+01, // 334
 6.835003560824457e-02, // 335
-3.524548959079112e+00, // 336
-1.152511713960376e+01, // 337
 2.497446390590568e-01, // 338
 1.677043950963147e+01, // 339
 1.154011354569148e-01, // 340
 7.414750123607935e+01, // 341
 3.653651712915493e-02, // 342
 4.224409165232077e+00, // 343
-4.020799518958716e-01, // 344
 9.571278665165149e-02, // 345
-1.453559856801310e+00, // 346
-1.437598747382498e+01, // 347
 2.014973381551670e+01, // 348
-4.437342158329273e+00, // 349
-3.714472529685672e-03, // 350
 9.083242214148776e-02, // 351
-1.996316209589908e+01, // 352
 4.100276694386930e-01, // 353
 1.811720061567939e+00, // 354
 1.828121308008920e+00, // 355
 3.696033610428785e+00, // 356
 2.216366907920183e-01, // 357
-4.867098613387440e+00, // 358
-1.157378368151389e+00, // 359
 3.537118986787274e-01, // 360
-2.218451185361990e-04, // 361
-5.579452759714679e+01, // 362
 1.695842937534920e-01, // 363
-5.145480959002310e-02, // 364
-5.703004577138220e-01, // 365
-2.863690782095930e-01, // 366
 1.093788292948795e+02, // 367
 9.235082784720072e+00, // 368
-2.367954144744845e+00, // 369
-1.774409270497165e-02, // 370
-1.913165725513419e+01, // 371
-4.612543554422313e-02, // 372
 6.756771427431080e+00, // 373
-5.311908166107838e+00, // 374
-2.674219390806348e+02, // 375
-8.499466513872513e+00, // 376
 3.042080826466585e+01, // 377
-4.302934491067051e+00, // 378
-2.963693990808501e-01, // 379
 1.306867462431039e+00, // 380
 7.608012363667174e+00, // 381
-1.794643790610261e+01, // 382
 2.376637839640942e+01, // 383
 6.507093844569205e+00, // 384
-1.162197992809380e+02, // 385
-9.432113017643460e-03, // 386
 4.892650801434193e+00, // 387
-1.823194610081394e-01, // 388
 7.618692459207208e+00, // 389
 3.700172196446437e-01, // 390
-2.055538811771168e+01, // 391
 4.338669474900285e-01, // 392
-3.937451589846184e-01, // 393
 3.622044615799072e+00, // 394
 3.065295309956497e+00, // 395
-8.135891695457099e-01, // 396
-3.994690089297591e+00, // 397
 5.654711340295648e+00, // 398
 8.652255664282996e-01, // 399
 1.005915265250811e+00, // 400
 1.823853590985876e+00, // 401
 1.687508520182778e+00, // 402
-3.299775433521350e+00, // 403
-1.044364877740423e-01, // 404
-9.908346188374371e-04, // 405
 3.058705120202319e-03, // 406
-6.611268295392998e+01, // 407
 2.842926970933198e-03, // 408
-1.386541714839225e+00, // 409
 3.567916629913991e+00, // 410
 3.933401413510115e-01, // 411
 6.740506754001553e-02, // 412
-1.947306572795300e+02, // 413
-8.783899796415757e-01, // 414
-2.855749043719053e-01, // 415
 2.669654786958036e+01, // 416
 1.174188232805364e+01, // 417
 5.846458808717615e+02, // 418
-3.096133674526233e+00, // 419
-6.749756955777836e+01, // 420
-4.027156399591018e+01, // 421
-6.019013261673378e-04, // 422
 4.025496380383903e+00, // 423
-9.954849442866216e-04, // 424
-1.118486767340564e-05, // 425
-1.741310986047283e-01, // 426
 1.569523831344603e+02, // 427
 2.264806113193226e+01 // 428
};

//----------------------------------------------------------------------------//

x2b_h2o_ion_v1x_p::x2b_h2o_ion_v1x_p()
{
    m_k_HH_intra =         2.151700564758605e-01; // A^(-1)
    m_k_OH_intra =         1.999999816503222e+00; // A^(-1)
                           
    m_k_XH_coul =          1.038128252838992e-01; // A^(-1)
    m_k_XO_coul =          9.196589412236915e-01; // A^(-1)
                           
    m_k_XLp_main =         5.960865361137950e-01; // A^(-1)
                           
    m_d_HH_intra =         1.894064210255164e+00; // A^(-1)
    m_d_OH_intra =         1.022971369863609e+00; // A^(-1)
                           
    m_d_XH_coul =          6.999993404082396e+00; // A^(-1)
    m_d_XO_coul =          5.399963541459254e+00; // A^(-1)
                           
    m_d_XLp_main =         6.310609579539030e+00; // A^(-1)
                           
    m_in_plane_gamma =     -9.721486914088159e-02;
    m_out_of_plane_gamma=  9.859272078406150e-02;

    m_r2i =  5.000000000000000e+00; // A
    m_r2f =  6.000000000000000e+00; // A

    std::copy(the_poly, the_poly + poly_type::size, m_poly);
};

//----------------------------------------------------------------------------//

static const x2b_h2o_ion_v1x_p the_model;

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

extern "C" {

//----------------------------------------------------------------------------//

#ifdef BGQ
void mbnrg_2b_h2o_li_poly(const double* w, const double* x ,
                    double* E, double* g1, double* g2)
#else
void mbnrg_2b_h2o_li_poly_(const double* w, const double* x ,
                    double* E, double* g1, double* g2)
#endif
{
    *E = the_model(w , x , g1, g2);
}

//----------------------------------------------------------------------------//

#ifdef BGQ
void mbnrg_2b_h2o_li_cutoff(double* r)
#else
void mbnrg_2b_h2o_li_cutoff_(double* r)
#endif
{
    *r = the_model.m_r2f;
}

//----------------------------------------------------------------------------//

} // extern "C"

////////////////////////////////////////////////////////////////////////////////
