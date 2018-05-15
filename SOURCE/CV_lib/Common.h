

/* define unit conversions */
#define BOHR_ANGS .529177
#define ANGS_BOHR 1.88972687777435527243
#define EV_HARTR 1.0/27.21138386
#define KV_HARTR 1.0/315777.0 
#define HARTR_KV 315777.0    /* Hartree to Kelvin */
#define TIME_CONV 0.0241888  /* dt/TIME_CONV; tau/TIME_CONV; femtosecond --> atomic unit */
#define PCONV 3.3989242e-09  /* Pext*PCONV;  atm --> atomic unit */
#define PROT_MASS 1822.8885  /* mass*PROT_MASS: atom mass --> atomic unit */


/* define constants */
#define MAX_Neighbour_Count 100
#define MaxAngMom 7


#define NINT(X) ((int) ((X)>=0.0 ? ((X)+0.5):((X)-0.5)))


typedef struct
{
  int Nsteps;       /* input before running */
  double dt;        /* timestep in PINY units */
  double Text;      /* External temp in Hartree */
  double THmat;     /* Hmatrix temp in Hartree */
  double Pext;      /* external pressure  */
  double tau_nhc;   /* NHC relaxation time */
  double tau_vol;   /* Hmatrix relaxation time */
  double tau_vol_nhc; /* Hmatrix NHC relaxation time */
  int NListUpdateStep; /* Neighborlist update timesteps */
}INPUT;


typedef struct
{
  double Re[3], Im[3];
} Global_CplxVec3;


typedef struct
{
  int nmax; /* max neighbors */
  int nlist[MAX_Neighbour_Count]; /* saves neighbor index */
  double ndist[3*MAX_Neighbour_Count]; /* saves dx, dy, dz in Angstroms */
  Global_CplxVec3 QlmDerivatives[MaxAngMom*MAX_Neighbour_Count];
  //  double *QlmDerivatives_Re, QlmDerivatives_Im*;
} NEIGHBORS;

typedef struct
{
  int *AtomIdx, AtomIdxMax;
  int CellNeighborIdx[27];
} Global_Cell;



typedef struct
{
  int CuDenLen, CuEmbedLen, CuPairLen;   /* EAM potential file lengths */
  char *CuDen, *CuPair, *CuEmbed;        /* EAM potential file names */
  double *DenSeparation, *ElecDensity, *ElecDenEmbed;  
  double *EnergyDenEmbed, *PairEnergy, *PairSeparation;

  gsl_interp_accel *acc1, *acc2, *acc3;
  gsl_spline *splineDen, *splineEmbed, *splinePair;
} POTSPLINE;


typedef struct
{
  int num_nhc;                 /* Num: # of NHC's  */
  int len_nhc;                 /* Num: Length of NHC's */
  int nres_nhc;                /* Num: # of RESPA NHC steps */
  int nyosh_nhc;               /* Num: # of Yosh NHC steps  */

  double dt_nhc,dti_nhc;       /* Num: NHC respa time steps  */
  double *wdti,*wdti2,*wdti4,*wdti8,*wdti16;   /* Lst: Yosh steps;  Lth:9 */

  double *mass_nhc,*gkt;       /* Lst: Mass,degs free of NHC's;  */
                               /* Lth: num_nhc x len_nhc  */

  double *x_nhc,*v_nhc,*f_nhc; /* Lst: Atm NHC pos,vel,for; */        
			       /* Lth: num_nhc x len_nhc    */
  double x_nhc_tot;
} THERM_INFO;


typedef struct
{
  double vol;
  double vol0;                 /* initial volume */

  int hmat_cons_typ;           /* Opt: Constraint on the box/h_mat */
			       /* 0=none,1=orthorhombic,2=monoclinic*/
  int iperd;

  double *hmato, *hmat_t, *hmat, *hmati;  /* Lst: Cell shape; Lth:9 */
  double mass_hm, text_hm;           /* Mass,  temperature of hh^-1 matrix */
  double c1_hm;                      /* Num: Useful constant */
  double *vgmat, *fgmat_p, *fgmat_v; /* Lst: Velocity, Force on hh^-1 matrix, Lth:9 */

  double *vtemps, *veigv,*veig;      /* Lst: Integration temp, Lth:9 */
  double *vexpdt, *vsindt;           /* Lst: Integration temp, Lth:9 */
  double *vtempx, *vtempv, *vtempf;  /* Lst: Integration temp, Lth:9 */

  double **x_vol_nhc, **v_vol_nhc;   /* Lst: Volume NHCs, Lth:len_nhc */
  double **f_vol_nhc, **mass_vol_nhc, **gkt_vol;  /* for massive nhc to h matrix */

  double *roll_mtx, *roll_mtv;
  double *roll_mtvv, *roll_mtf;
} PAR_RAHMAN;


typedef struct
{
  double tvten[9];                  /* Lst: KE  tensor           ; Lth: 9  */
  double pvten_Global[9];           /* Lst: PV  tensors          ; Lth: 9  */
  double pvten_tmp[9];              /* Lst: PV  tensors for SHAKE; Lth: 9  */
} PTENSOR;


typedef struct
{
  double kappa_ev[2];               /* Lst : coupling k        */
  double x_ev[2];                   /* Lst : pos of EVs        */
  double f_ev[2];                   /* Lst : force of EVs      */  
} EXT_VAR;


typedef struct
{
  int ind_ev;                       /* index of the extended variable */
  double x_cv;                      /* pos. of collective var.     */
} COLLC_VAR;


typedef struct
{
  char *symb;
  double *pos, *vel, *mass, *force;
}CLATOMS;


typedef struct
{
  int l_ordr;                       /* order of harm. spher. */

  double depth_nb;                  /* skin depth of the neighborhood list */
  double rmin_ql, rmax_ql;          /* upper & lower cutoff for switch function */

  double *cos_mphi, *sin_mphi;      /* Lth: l_ordr+1, cos(m*phi), sin(m*phi) and 0<=m<=l_ordr */
  double *dcos_mphi_x, *dsin_mphi_x;/* Lth: l_ordr+1, d(cos(m*phi))/dx1, d(sin(m*phi))/dx1, 0<=m<=l_ordr */
  double *dcos_mphi_y, *dsin_mphi_y;/* Lth: l_ordr+1, d(cos(m*phi))/dy1, d(sin(m*phi))/dy1, 0<=m<=l_ordr */

  double *alp;                      /* Lth: l_ordr+1, associate legendre polynomial(ALP) */
  double *alp_d;                    /* Lth: l_ordr+1, first derivative of (ALP)  */
  double *qlm_re, *qlm_im;          /* Lth: l_ordr+1, Qlm real, imaginary */
  double *pre_sph2;                 /* Lth: l_ordr+1, prefactor of sph harm */

  double *sh_re, *sh_im;            /* Lth: l_ordr+1  real, imaginary of sph harm */
} STEINHDT_PKG;


typedef struct
{
  double kinet, vpot_emd;             /* kinetic energy, potential energy of particles */
  double vpotnhc, kinet_nhc;      /* potential and kinetic energy of nhc of particles */
  double kinet_nhc_v, vpotnhc_v;  /* kinetic and potential energy of nhc of h matix */ 
  double kinet_v, vpot_v;         /* kinetic and potential energy for h matrix */    
  double vpot_Q6,vpot_Q4;
  double vpot;
} STAT_AVG;


