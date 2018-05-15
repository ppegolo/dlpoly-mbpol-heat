
double Steinhardt_Ql(int natoms, int l_ordr, double *ev_pos, 
		     double *ev_force, double *ev_kappa, double *cv_pos, 
		     double RMIN_Ql, double RMAX_Ql, NEIGHBORS *NList, 
		     int cv_idx, double *force, double *pvtens)
{
  /* Steinhardt orientation order parameter calculation using GSL */


  /*
    1. pos: pointer to an array (3*Natoms)
    2. Natoms: number of atoms
    3. l_ordr: 4 or 6, corresponds to Q4 or Q6
    4. ev_pos: position of extended variables, array
    5. ev_force: forces on the extended variables, array
    6. ev_kappa: coupling constants, array
    7. cv_pos: position of collective variables
    8. RMAX_Ql, RMIN_Ql: min and max cutoff values for Ql calculation
    9. NEIGHBORS: array of structures
    10. cv_idx: index of collective/extended variable
    11. force: 3*Natoms, forces on the atoms/ions
    12. pvten: 9x1 vector, pressure tensor
   */

  int i, j, k, Ni;
  
  double arg, arg1, f_pre, comm_fac, DX[3], dr, idr, idr3;
  double r12xy, r12xy_inv, r12xy_inv3, dx[3];
  double cos_th, dcos_th_x, dcos_th_y, dcos_th_z, fc, dfc[3], ql_v;
  double tmpalp[l_ordr+1], tmpalp_d[l_ordr+1], ForceLoc[3], diff, Energy_ExtVar;
  double pvten_tmp[9];
  
  double *pre_sph2 = NULL, *alp = NULL, *alp_d = NULL;
  double *cos_mphi = NULL, *sin_mphi = NULL;

  double *dcos_mphi_x = NULL, *dsin_mphi_x = NULL;
  double *dcos_mphi_y = NULL, *dsin_mphi_y = NULL;

  double *qlm_re = NULL, *qlm_im = NULL;
  double *sh_re = NULL, *sh_im = NULL;


  /* START : Memory allocation */
  qlm_re = (double *) malloc ((l_ordr+1)*sizeof(double));
  qlm_im = (double *) malloc ((l_ordr+1)*sizeof(double));

  pre_sph2 = (double *) malloc ((l_ordr+1)*sizeof(double));

  sh_re = (double *) malloc ((l_ordr+1)*sizeof(double));
  sh_im = (double *) malloc ((l_ordr+1)*sizeof(double));

  cos_mphi = (double *) malloc ((l_ordr+1)*sizeof(double));
  sin_mphi = (double *) malloc ((l_ordr+1)*sizeof(double));

  dcos_mphi_x = (double *) malloc ((l_ordr+1)*sizeof(double));
  dsin_mphi_x = (double *) malloc ((l_ordr+1)*sizeof(double));
  dcos_mphi_y = (double *) malloc ((l_ordr+1)*sizeof(double));
  dsin_mphi_y = (double *) malloc ((l_ordr+1)*sizeof(double));

  alp = (double *) malloc ((l_ordr+1)*sizeof(double));
  alp_d = (double *) malloc ((l_ordr+1)*sizeof(double));
  /* END : Memory allocation */


  /* START : Initialize variables */
  for(i = 0; i <= l_ordr; i++)
    {
      qlm_re[i] = 0.00; 
      qlm_im[i] = 0.00; 
    } 

  pre_sph2[0] = 1.0*0.5;  /* divide by 2, to be required later */
  for(i = 1; i <= l_ordr; i++)
    {
      pre_sph2[i] = 1.0;
      for(j = 1; j <= 2*i; j++)
	pre_sph2[i] *= (l_ordr-i+j);

      pre_sph2[i] = 1.0/pre_sph2[i];
    }
  /* END : Initialize variables */  


  /* START : Ql calculation */  
  for( i = 0; i <= (natoms-2); i++)
    { /* the last atom has an empty neighbor list */
      for (Ni = 0; Ni < NList[i].nmax; Ni++)
        {
          j = NList[i].nlist[Ni];
	  
	  /* dx is in Angstroms */
	  dx[0] = -1*(NList[i].ndist[3*Ni+0]);
	  dx[1] = -1*(NList[i].ndist[3*Ni+1]);
	  dx[2] = -1*(NList[i].ndist[3*Ni+2]);
	  
	  dr = NormV3(dx);  /* dr in Angstroms */
	  idr = 1.0/dr;
	  idr3 = idr*idr*idr;

	  r12xy = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
	  r12xy_inv = 1.0/r12xy;
	  r12xy_inv3 = r12xy_inv*r12xy_inv*r12xy_inv;

	  /*  smoothing function */
	  if(dr <= RMIN_Ql)
	    {
	      fc = 1.0;
	      InitV3(dfc);
	    }
	  else if(dr > RMAX_Ql)
	    {
	      fc = 0.0;
	      InitV3(dfc);
	    }
	  else
	    {
	      arg1 = M_PI/(RMAX_Ql - RMIN_Ql);
	      fc = 0.5*cos((dr - RMIN_Ql)*arg1) + 0.5; /* smoothing function */

	      /* derivative of the smoothing function */
	      arg = -0.5*arg1*sin((dr - RMIN_Ql)*arg1)*idr;
	      V3SmulV3(dx, -arg, dfc);
	    }

	  /* Ql calculation */
	  if(dr < RMAX_Ql)
	    {
	      cos_th    = dx[2]*idr;
	      
	      /* cos(mPhi), m = 0, 1 */
	      cos_mphi[0] = 1.0;
	      cos_mphi[1] = dx[0]*r12xy_inv;

	      /* sin(mPhi), m = 0, 1 */
	      sin_mphi[0] = 0.0;
	      sin_mphi[1] = dx[1]*r12xy_inv;

	      /* derivative cos(mPhi), m = 0, 1 */
	      dcos_mphi_x[0] = 0.0;
	      dcos_mphi_x[1] = dx[0]*dx[0]*r12xy_inv3 - r12xy_inv;
	      dcos_mphi_y[0] = 0.0;
	      dcos_mphi_y[1] = dx[0]*dx[1]*r12xy_inv3;

	      /* derivative sin(mPhi), m = 0, 1 */
	      dsin_mphi_x[0] = 0.0;
	      dsin_mphi_x[1] = dcos_mphi_y[1];
	      dsin_mphi_y[0] = 0.0;
	      dsin_mphi_y[1] = dx[1]*dx[1]*r12xy_inv3 - r12xy_inv;

	      /* obtain cos(mPhi) and sin(mPhi) for m >= 2 */
	      for(k = 2; k <= l_ordr; k++)
		{
		  cos_mphi[k] = 2.0*cos_mphi[1]*cos_mphi[k-1] - cos_mphi[k-2];
		  sin_mphi[k] = 2.0*cos_mphi[1]*sin_mphi[k-1] - sin_mphi[k-2];

		  /* derivatives */
		  dcos_mphi_x[k]  = 2.0*(dcos_mphi_x[1]*cos_mphi[k-1]
					 +cos_mphi[1]*dcos_mphi_x[k-1]) - dcos_mphi_x[k-2];
		  dcos_mphi_y[k]  = 2.0*(dcos_mphi_y[1]*cos_mphi[k-1]
					 +cos_mphi[1]*dcos_mphi_y[k-1]) - dcos_mphi_y[k-2];
		  dsin_mphi_x[k]  = 2.0*(dcos_mphi_x[1]*sin_mphi[k-1]
					 +cos_mphi[1]*dsin_mphi_x[k-1]) - dsin_mphi_x[k-2];
		  dsin_mphi_y[k]  = 2.0*(dcos_mphi_y[1]*sin_mphi[k-1]
					 +cos_mphi[1]*dsin_mphi_y[k-1]) - dsin_mphi_y[k-2];
		}
	      
	      /* use GSL to calculate the correct associated Legendre polynomials */
	      for(k = 0; k <= l_ordr; k++)
		{
		  gsl_sf_legendre_Plm_deriv_array(l_ordr, k, cos_th, tmpalp, tmpalp_d);
		  alp[k] = tmpalp[l_ordr-k];
		  alp_d[k] = tmpalp_d[l_ordr-k];
		}

	      for(k = 0; k <= l_ordr; k++)
		{
		  dcos_th_x = dx[2]*dx[0]*idr3;
		  dcos_th_y = dx[2]*dx[1]*idr3;
		  dcos_th_z = dx[2]*dx[2]*idr3 - idr;
		  
		  /* spherical Harmionics Real, Im */
		  sh_re[k] = cos_mphi[k]*alp[k];
		  sh_im[k] = sin_mphi[k]*alp[k];
		  
		  /* ql real, imaginary */
		  qlm_re[k] += fc*sh_re[k];  
		  qlm_im[k] += fc*sh_im[k];

		  /* derivatives w.r.t particle positions */
		  NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[0] = sh_re[k]*dfc[0] + 
		    fc*(cos_mphi[k]*dcos_th_x*alp_d[k]+alp[k]*dcos_mphi_x[k]);
		  NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[1] = sh_re[k]*dfc[1] + 
		    fc*(cos_mphi[k]*dcos_th_y*alp_d[k]+alp[k]*dcos_mphi_y[k]);
		  NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[2] = sh_re[k]*dfc[2] + 
		    fc*(cos_mphi[k]*dcos_th_z*alp_d[k]);

		  NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[0] = sh_im[k]*dfc[0] + 
		    fc*(sin_mphi[k]*dcos_th_x*alp_d[k]+alp[k]*dsin_mphi_x[k]);
		  NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[1] = sh_im[k]*dfc[1] + 
		    fc*(sin_mphi[k]*dcos_th_y*alp_d[k]+alp[k]*dsin_mphi_y[k]);
		  NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[2] = sh_im[k]*dfc[2] + 
		    fc*(sin_mphi[k]*dcos_th_z*alp_d[k]);
		}
	    } /* endif: exclude pair distance > rmax_ql*/
	} /* endfor Ni loop */
    } /* endfor i loop */
  /* END : Ql calculation */

  comm_fac  = 1.0/((double) (natoms*natoms));

  ql_v = 0.0;
  for(i = 0; i <= l_ordr; i++)
  {  ql_v += pre_sph2[i]*(qlm_re[i]*qlm_re[i] + qlm_im[i]*qlm_im[i]);
    //     printf("qlm %d %.15g %.15g\n",i,qlm_re[i],qlm_im[i]);
  }
  ql_v = sqrt(ql_v*2.0*comm_fac);  // using symmetry: |Qlm|^2 = |Ql,-m|^2
  
  cv_pos[cv_idx] = ql_v;
  //  printf("The Q%d for the Cfg file is %f\n", l_ordr, ql_v);

  /* potential/force for extended variable */
  diff = ev_pos[cv_idx] - ql_v;   /* (ExtVar6-Q6), (ExtVar4-Q4), are dimensionless */
  ev_force[cv_idx]  = -ev_kappa[cv_idx]*diff;   /* force in Hartree */
  Energy_ExtVar = ev_kappa[cv_idx]*diff*diff*0.5;   /* energy in same units as kappa */
  //  printf("The force and energy on ExtVar%d are %f, %f\n", l_ordr, ext_var.f_ev[Ql->ind_ev], Energy_ExtVar);

  /*=======================================================================*/
  /*  iii) get forces of each pair */
  f_pre = comm_fac/ql_v;

  InitV9(pvten_tmp);
  
  for( i = 0; i <= (natoms-2); i++)
    {     /* the last atom has an empty neighbor list */
      for (Ni = 0; Ni < NList[i].nmax; Ni++)
        {
          j = NList[i].nlist[Ni];

	  DX[0] = -1*(NList[i].ndist[3*Ni+0]);
	  DX[1] = -1*(NList[i].ndist[3*Ni+1]);
	  DX[2] = -1*(NList[i].ndist[3*Ni+2]);

          dr = (NormV3(DX));

	  ForceLoc[0] = 0.00;
	  ForceLoc[1] = 0.00;
	  ForceLoc[2] = 0.00;
	  
	  if(dr < RMAX_Ql)  /* dr, rmax_ql in Bohr */
	    {
	      for(k = 0; k <= l_ordr; k++)
		{
		  ForceLoc[0] += pre_sph2[k]*(qlm_re[k]*NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[0]
					      + qlm_im[k]*NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[0]);

		  ForceLoc[1] += pre_sph2[k]*(qlm_re[k]*NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[1]
					      + qlm_im[k]*NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[1]);

		  ForceLoc[2] += pre_sph2[k]*(qlm_re[k]*NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Re[2]
					      + qlm_im[k]*NList[i].QlmDerivatives[(l_ordr+1)*Ni+k].Im[2]);
		}  /* endfor k loop */
	      /* times 2 because of symmetry: ReQlm [ReQlm]' = ReQl,-m [ReQl,-m]' */
	      
	      ForceLoc[0] *= 2.0*f_pre*(ev_force[cv_idx]);
	      ForceLoc[1] *= 2.0*f_pre*(ev_force[cv_idx]);
	      ForceLoc[2] *= 2.0*f_pre*(ev_force[cv_idx]);

	    } /* endif */

	  force[3*i+0] += ForceLoc[0];
	  force[3*i+1] += ForceLoc[1];
	  force[3*i+2] += ForceLoc[2];

	  force[3*j+0] += -ForceLoc[0];
	  force[3*j+1] += -ForceLoc[1];
	  force[3*j+2] += -ForceLoc[2];
	  
          DX[0] = -1*DX[0];
          DX[1] = -1*DX[1];
          DX[2] = -1*DX[2];
	  V3V3diadicV9addV9(DX, ForceLoc, pvten_tmp);

	} /* end for Ni */
    } /* end for i */

  /* pressure tensor calculations */
  V9symmetrizeV9addV9(pvten_tmp, pvtens);

  free(qlm_re);
  free(qlm_im);
  free(pre_sph2);
  free(sh_re);
  free(sh_im);
  free(cos_mphi);
  free(sin_mphi);
  free(dcos_mphi_x);
  free(dsin_mphi_x);
  free(dcos_mphi_y);
  free(dsin_mphi_y);
  free(alp);
  free(alp_d);

  return Energy_ExtVar;
}
