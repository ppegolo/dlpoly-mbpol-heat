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

//fit-fullpolargrid-fixedwaterparams 100%polfac 03/15/17
static const double the_poly[] = {
 1.118689224011763e+02, // 0
-1.558086633268414e+03, // 1
 1.092189792217720e+03, // 2
 1.111234488756074e+01, // 3
 7.270784361491772e+00, // 4
 2.378907402563357e+02, // 5
 1.855001355937474e+03, // 6
-6.874536648965535e+00, // 7
 1.228118171562922e+02, // 8
-1.285300792925100e+03, // 9
-3.223368028603641e+01, // 10
-3.531280907153072e+00, // 11
-3.303559119396961e+02, // 12
 6.776240865888720e+00, // 13
-2.807544218173088e+02, // 14
 6.678723538889631e+02, // 15
-4.034619235285903e+02, // 16
 2.727823966324587e+01, // 17
 1.446418976066514e+01, // 18
 3.551117952873494e+02, // 19
 1.537335988673700e+02, // 20
 5.821441307513138e+01, // 21
 1.378799162602521e+01, // 22
 7.127094787985552e+00, // 23
-6.781797297919020e+01, // 24
 1.358449315491242e+01, // 25
 2.161018921497626e+02, // 26
 9.980744451524498e+02, // 27
-3.595815212265356e+02, // 28
-3.278378179010325e+02, // 29
 3.975713057078425e+00, // 30
 3.524010569787907e+01, // 31
 1.242038759123033e+00, // 32
 1.373786504559208e+02, // 33
-1.530437509933577e+03, // 34
 4.071211460676246e+02, // 35
-7.272232753070238e-01, // 36
-6.518090731724902e+01, // 37
-5.737596678026878e+02, // 38
-4.872935719142394e+02, // 39
-4.041821071056260e+02, // 40
-3.335854939147913e+00, // 41
-1.141996761086143e+01, // 42
 5.647967540037731e-02, // 43
-3.647251363712593e+02, // 44
-3.832812213705397e+01, // 45
 3.324247823667760e+01, // 46
 1.184912162415162e+01, // 47
-3.195633536092347e+01, // 48
 3.285353859606557e+00, // 49
-4.268824892788059e+01, // 50
-2.115235139536318e+01, // 51
 3.078316381648446e+01, // 52
 4.535361719111894e+01, // 53
 6.038553365481836e+02, // 54
 3.699126423959425e+02, // 55
 3.960563431205756e+01, // 56
 5.870206222464907e+02, // 57
-1.236341985788662e+02, // 58
-1.856377750501828e+02, // 59
-2.097175074825028e+00, // 60
-1.577043010608390e+02, // 61
 1.902219057647398e-02, // 62
-1.773265060149170e+00, // 63
-8.969867954780071e+00, // 64
-1.980139801828772e+01, // 65
-7.063261180249685e+01, // 66
-4.137351344313408e+00, // 67
 2.889586672588148e+02, // 68
-4.571807169198219e-04, // 69
-5.498150057216998e+01, // 70
 2.636471574276082e-02, // 71
 1.832830115655273e+01, // 72
-9.598883870937774e+02, // 73
 1.147741694288335e+00, // 74
-1.832671263467880e+02, // 75
-1.355005176825786e+01, // 76
 6.456032867531039e+00, // 77
-9.962700461609088e-01, // 78
 4.227657869714574e+01, // 79
 2.270074172187727e-03, // 80
-3.853435169682438e+02, // 81
-5.371500181856678e-03, // 82
 1.803268214279502e-03, // 83
-2.655792812140274e-03, // 84
-1.641432738850801e-05, // 85
 1.632306412641268e+01, // 86
 9.295675373547588e-01, // 87
 1.136645658600747e+03, // 88
 3.501840570906024e-01, // 89
-3.835104350570389e+00, // 90
 1.282811891069206e+02, // 91
-1.002231768145675e+01, // 92
-8.993600785415319e+01, // 93
 3.507674063489011e+01, // 94
 2.217880807522698e+02, // 95
-1.800721520124600e+00, // 96
 2.408937436759236e+01, // 97
-3.923728324371851e+02, // 98
-2.825800923597247e+00, // 99
-2.567129602490528e-01, // 100
 1.416948117668975e-03, // 101
-8.430515908172025e-01, // 102
-9.643030547479788e+01, // 103
-2.733978868105109e+01, // 104
-3.154767116255481e+00, // 105
-1.488816620580689e-01, // 106
-1.819991797729027e+02, // 107
-9.015714332797918e+01, // 108
-2.396547923211739e+00, // 109
-4.717704498863527e-01, // 110
 4.640044985711294e+00, // 111
-6.019193786384934e+00, // 112
 4.145785239026792e-04, // 113
-2.905518709962152e-02, // 114
-3.252676209540530e+02, // 115
-7.242039289238767e+01, // 116
 7.641739328898355e+01, // 117
-1.094295327669893e+00, // 118
-7.426950952055287e+00, // 119
 1.503707951451348e+01, // 120
 1.165849888037632e+02, // 121
 5.519269111549092e-01, // 122
 5.789318280254636e+01, // 123
 6.192260967252062e+02, // 124
-1.794722322725870e+00, // 125
 2.983945141802863e-01, // 126
-3.155686174620947e+02, // 127
 7.411374799100902e+01, // 128
-1.891189537355061e-02, // 129
 1.097960572725427e+02, // 130
-5.373253567395082e+02, // 131
-1.176416600852098e+01, // 132
-5.632207972310768e+00, // 133
 2.117931267985845e+01, // 134
 1.524089202793523e+00, // 135
 3.595203193794744e-02, // 136
-5.791438928320870e-02, // 137
-3.931360423179168e+01, // 138
 6.242213446054161e+02, // 139
-1.449852246241461e-02, // 140
-1.328958038347467e+01, // 141
 1.108176935972844e+02, // 142
-4.738469365754408e+02, // 143
 1.781468342001212e+01, // 144
-3.805019477145925e-01, // 145
-2.719396993589112e+01, // 146
-1.792574448896592e+00, // 147
 2.824369317820438e+02, // 148
 4.669130181710096e+01, // 149
-2.892390914939843e+01, // 150
 4.870244704117164e+00, // 151
 1.977825775929892e+02, // 152
 2.403310404127136e-01, // 153
-1.117425749660582e-01, // 154
 4.948225408128215e-01, // 155
 6.059184893828735e+02, // 156
 2.381152804170414e+00, // 157
-4.760233539261344e-02, // 158
-2.189689508876273e-02, // 159
 5.923077242623447e+00, // 160
-3.796997129774886e-05, // 161
-4.416966485364160e-01, // 162
 1.109752471796873e+01, // 163
 2.141634754483150e+01, // 164
 9.977291937521078e-01, // 165
-5.811095550740085e+01, // 166
-3.508294767478670e+00, // 167
 4.649330979204934e-01, // 168
-1.208508189765898e+00, // 169
-4.240338595271390e+01, // 170
 1.782561497567371e-03, // 171
 5.669679777149897e+00, // 172
 1.812425808151612e+00, // 173
-1.017137187035662e+01, // 174
-1.850495187326429e+02, // 175
-6.280230074347437e-01, // 176
 2.784545063443423e-02, // 177
 8.556785120243134e+00, // 178
 1.085550275580849e-01, // 179
 1.431731897725528e+01, // 180
-4.444644361787099e+02, // 181
 3.307024608090718e+00, // 182
-1.949314962943324e-02, // 183
 7.524415712307886e+02, // 184
 3.337047966385399e-05, // 185
 2.349611913673884e-04, // 186
 4.551894810356409e+01, // 187
-4.296475233512196e-04, // 188
 1.135473343673134e-01, // 189
-1.043834438811577e+01, // 190
-8.133744392836685e+02, // 191
-1.921062381663353e-01, // 192
 6.404136869901346e-01, // 193
-3.436263302755802e+00, // 194
-1.540620152937419e-03, // 195
 8.377372602005304e-02, // 196
 4.866748625972376e-03, // 197
 6.886258371531514e+01, // 198
 3.846118977484519e-01, // 199
-8.597261634405430e-01, // 200
-7.711587708870794e-03, // 201
 1.338277133604781e-01, // 202
 1.348151823593798e+01, // 203
 2.781772640360290e-01, // 204
-1.847032433821116e+00, // 205
-3.342517772550885e+01, // 206
 6.637088154702095e-04, // 207
-2.224667114941217e+02, // 208
-2.130101541489825e+01, // 209
-3.181633514864437e-04, // 210
 6.679786656732380e-01, // 211
-1.106026623011486e-06, // 212
-4.950750208098255e-02, // 213
-9.588384291360303e+00, // 214
-6.012368751942745e-02, // 215
 1.660275356658160e+00, // 216
-7.730754254836804e-01, // 217
 6.014585880316287e+00, // 218
 8.046656883628544e-04, // 219
 7.264227655657140e-04, // 220
 2.460380130867961e-01, // 221
-1.734776758089639e+01, // 222
 3.199644034443905e+00, // 223
-5.774551648288261e-03, // 224
-5.520454954716038e+01, // 225
-1.750177975638350e-04, // 226
-2.388580258410326e-08, // 227
-3.281810021742582e-03, // 228
 1.528543290795386e-02, // 229
 3.677868586326349e-01, // 230
-7.252029253866573e+00, // 231
 9.423090830031392e-01, // 232
-3.647151816635629e-05, // 233
 5.060940860713563e+01, // 234
-4.155242533718199e-02, // 235
 1.469585302596687e+01, // 236
-2.460333484996169e+01, // 237
 1.715480485479799e-01, // 238
-4.477697894358992e+01, // 239
-1.156726725441348e+01, // 240
 9.087413668811541e-02, // 241
-6.768856034060270e-03, // 242
 2.162215934296337e+00, // 243
-3.723695075192654e-06, // 244
 1.800683270430989e-03, // 245
 1.361937910078172e+01, // 246
-2.142698950367774e-04, // 247
 1.681363114976093e-03, // 248
-3.703117094471524e-03, // 249
-2.621620392908395e-04, // 250
-6.480681967365171e-01, // 251
 1.849131816792689e-01, // 252
-2.315381757241610e-02, // 253
 5.833083891461976e+01, // 254
 3.337299113293032e-01, // 255
 9.284143697871654e-03, // 256
 2.378899203163833e+01, // 257
-4.558842288897573e-02, // 258
 2.524631073167432e+02, // 259
 2.261447438384827e+02, // 260
 1.300238000478088e+00, // 261
-1.786350271804076e+01, // 262
 3.122426126211358e-01, // 263
-8.095666498913049e-01, // 264
-6.703197550174687e-02, // 265
 4.731761333054564e-04, // 266
-5.373625713848036e-01, // 267
-1.521700192831657e+00, // 268
-1.157205136321357e-01, // 269
 4.064640783384400e-02, // 270
-1.535004578138493e+00, // 271
 1.941934563171621e+02, // 272
-1.588299271698038e+01, // 273
-3.501452386130351e+00, // 274
-5.269227769129647e-02, // 275
-1.854628371579154e+00, // 276
-1.160527495174238e-03, // 277
-7.854346077468454e-01, // 278
 7.474831507026903e-03, // 279
 1.798571928216067e-05, // 280
-1.620420204979605e+01, // 281
-1.530531255586894e+00, // 282
-7.301836404813048e+00, // 283
-2.819535503090102e-03, // 284
 5.222470600762092e-01, // 285
-9.622117516839469e-02, // 286
-3.241944761374119e-02, // 287
 7.428434013939099e+00, // 288
-9.721270974888728e-06, // 289
 4.226321269619336e-04, // 290
 3.035100279083159e+01, // 291
-3.650286487765725e+01, // 292
 1.777181010299608e+00, // 293
-8.546243600429939e+01, // 294
 7.904148470829188e-01, // 295
-1.935161012823615e-01, // 296
 7.923108340184056e+00, // 297
-4.593381257369033e-03, // 298
 3.941150592226066e+00, // 299
 2.069586210696164e+00, // 300
-1.734660013776785e+00, // 301
 1.038795553107363e-07, // 302
-2.050671175610478e-01, // 303
 1.058996323855230e-03, // 304
 1.390521533826182e+00, // 305
 4.050171250853716e+01, // 306
-6.957530923464111e-07, // 307
 2.108589331462328e+00, // 308
 5.729011655171536e-08, // 309
-6.838393535154663e+01, // 310
-1.139000399210430e+00, // 311
-1.294475839479995e+01, // 312
 2.489601456430276e+01, // 313
-6.130904322377356e-01, // 314
-7.734531816207568e-01, // 315
 1.359482751448900e-06, // 316
 2.244482096238645e+01, // 317
-3.237085993866357e-06, // 318
 1.793637306562432e-03, // 319
 2.896922894171731e+01, // 320
 1.984849518891021e+01, // 321
 1.615745384572521e-03, // 322
-5.369163426040357e-02, // 323
-1.734270625480452e+02, // 324
 5.833564185338579e+00, // 325
-1.632482838506342e-01, // 326
 4.585684552137376e+01, // 327
 2.799705170973102e+00, // 328
 6.743061622812293e-03, // 329
 1.113818238514249e-02, // 330
 1.521772622666250e+01, // 331
-9.158791463111970e-02, // 332
 3.966973143110358e+00, // 333
 2.175894099786609e-01, // 334
-2.678842761756216e-07, // 335
-5.701953270368211e-03, // 336
 2.969188635824135e-06, // 337
-4.238068482642115e+01, // 338
-5.979923061941063e+00, // 339
-1.051904482662608e+01, // 340
-9.568850755377399e-04, // 341
-4.377076315844065e-01, // 342
 4.727557874568487e+00, // 343
 3.649879349413231e-02, // 344
 4.785791757555918e+01, // 345
 2.667649176906808e+00, // 346
-1.259417585961010e-03, // 347
 9.940947930289261e+01, // 348
 4.102004182435931e-04, // 349
 3.564732219872449e+00, // 350
-7.923839956002692e-01, // 351
 1.861059469427351e+00, // 352
-6.142069536614079e-03, // 353
 7.101045695852501e-06, // 354
 4.667091403042205e+02, // 355
 1.256719805541106e+00, // 356
 1.360572815001792e-05, // 357
 1.109286891670170e+00, // 358
-1.147485656765055e-01, // 359
 2.944610393790537e-01, // 360
 2.494464877456725e-01, // 361
 6.720808713486982e-04, // 362
-4.645357565020824e+01, // 363
 3.180092395018285e+00, // 364
-2.937964341114079e+01, // 365
-2.838353551878574e-03, // 366
-9.345442843079949e-03, // 367
-1.415164607936164e-01, // 368
 2.608867207957114e+00, // 369
 7.044765040459624e-05, // 370
 7.401890377977845e+02, // 371
 4.288768079177214e-02, // 372
-2.765822264924603e-05, // 373
 1.666119681618496e-01, // 374
 3.716506277072302e-01, // 375
-4.790484556187814e-01, // 376
 2.442440394046809e-03, // 377
-2.787264358601726e+01, // 378
 4.978246594767569e+01, // 379
-5.689150426339891e-04, // 380
-9.291069456549826e+00, // 381
 2.227167208630235e-05, // 382
 2.821597411960386e-05, // 383
-1.319954133057158e+03, // 384
 8.336693129256864e-04, // 385
 3.125044721544301e-04, // 386
-3.446921067650163e+02, // 387
 4.932011097525909e+00, // 388
-6.384843231197498e-02, // 389
-1.473459106244696e+02, // 390
 1.984272200387490e-02, // 391
-1.589779599332500e-01, // 392
 1.298317938403327e+01, // 393
-1.334070502723418e-05, // 394
 5.068533989323892e-03, // 395
-1.403841780142196e+00, // 396
 2.503167867790113e-02, // 397
 7.013536975404150e-01, // 398
-3.735009729686770e+00, // 399
 9.478516258063903e-01, // 400
-3.638735800711775e+01, // 401
 1.026818797133086e-06, // 402
-5.220763293609670e-02, // 403
-8.310094983978347e-01, // 404
-4.373062832237252e+00, // 405
-7.029016493903015e+00, // 406
-8.999371975855763e-04, // 407
-3.388427769182683e+00, // 408
 2.817998250437343e-03, // 409
 3.860346571115308e+01, // 410
-1.462280627184749e+01, // 411
-1.884916604871195e-01, // 412
 1.878195255411979e+02, // 413
 2.854605628594311e+01, // 414
-6.742017594241582e-01, // 415
-6.441097991494068e-04, // 416
-5.519211585293593e-01, // 417
-2.482653184244645e-01, // 418
 3.500498245345825e-02, // 419
-5.421784089130426e-04, // 420
 1.325160863170309e-02, // 421
 3.669682239241531e-01, // 422
-5.789167270918158e-05, // 423
-6.026579873936947e-01, // 424
 1.534261042268642e-02, // 425
 2.118736645242532e+00, // 426
 5.618337217506856e-02, // 427
 7.350825134957270e-03  // 428
};

//----------------------------------------------------------------------------//

x2b_h2o_ion_v1x_p::x2b_h2o_ion_v1x_p()
{
    m_k_HH_intra =         2.101776211312097e-01; // A^(-1)
    m_k_OH_intra =         3.033739168595145e-01; // A^(-1)
                           
    m_k_XH_coul =          8.514971338016590e-01; // A^(-1)
    m_k_XO_coul =          8.914473425594043e-01; // A^(-1)
                           
    m_k_XLp_main =         1.039449866720506e+00; // A^(-1)
                           
    m_d_HH_intra =         6.366721884964642e-01; // A^(-1)
    m_d_OH_intra =         1.443106669594841e+00; // A^(-1)
                           
    m_d_XH_coul =          6.189777918308249e+00; // A^(-1)
    m_d_XO_coul =          6.238568568191817e+00; // A^(-1)
                           
    m_d_XLp_main =         3.909368547725790e+00; // A^(-1)
                           
    m_in_plane_gamma =     -9.721486914088159e-02;
     m_out_of_plane_gamma= 9.859272078406150e-02;

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
void mbnrg_2b_h2o_f_poly(const double* w, const double* x ,
                    double* E, double* g1, double* g2)
#else
void mbnrg_2b_h2o_f_poly_(const double* w, const double* x ,
                    double* E, double* g1, double* g2)
#endif
{
    *E = the_model(w , x , g1, g2);
}

//----------------------------------------------------------------------------//

#ifdef BGQ
void mbnrg_2b_h2o_f_cutoff(double* r)
#else
void mbnrg_2b_h2o_f_cutoff_(double* r)
#endif
{
    *r = the_model.m_r2f;
}

//----------------------------------------------------------------------------//

} // extern "C"

////////////////////////////////////////////////////////////////////////////////