#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_legendre.h>
#define max_ncount 100
#define max_angmom 7


typedef struct
{
  double Re[3], Im[3];
} Global_CplxVec3;

typedef struct
{
  int nmax; /* max neighbors */
  int nlist[max_ncount]; /* saves neighbor index */
  double ndist[3*max_ncount]; /* saves dx, dy, dz in Angstroms */
  Global_CplxVec3 QlmDerivatives[max_angmom*max_ncount];
} NEIGHBORS;

#include "Common.c"
#include "Nlist.c"
#include "Ql.c"

// very inefficient (but easy) for now, 
// will improve if this becomes slow enough to be noticeable
void  ql_wrapper_(int* comm_bead, int* bead_rank,int* bead_size,
                  int* np, double* NList_Cutoff, double* Hmat, double *posX, double *posY, double *posZ, 
		  double *ev_pos, double *ev_force, double *ev_kappa, double *cv_pos, double* rmin_ql, 
		  double* rmax_ql, double *Force, double *pvtens, double *Ql)
{
  double H[3][3], *pos = NULL, *force = NULL;
  int i, j, k, Natoms;

  // no parallelization for now
  if(*bead_rank != 0) return; 

  for(i = 0, k = 0; i < 3; i++) 
    for(j = 0; j < 3; j++, k++) 
      H[i][j] = Hmat[k];


  Natoms = (int) (*np/3.0); // number of oxygen atoms 

  force = (double *) malloc(3*Natoms*sizeof(double));
  pos =  (double *) malloc(3*Natoms*sizeof(double));
  //  pvtens = (double *) malloc(9*sizeof(double));

  //    IMPORTANT: the atomic positions are assumed to be in the order: O H H O H H O H H ... 
  for (i = 0; i < Natoms; i++)
    {
      j = 3*i;

      // pos saves positions of ONLY oxygen atoms
      pos[3*i+0] = posX[j];   /* posX, pos in cartesian coordinates, Angstroms */
      pos[3*i+1] = posY[j];   /* posY, pos in cartesian coordinates, Angstroms */
      pos[3*i+2] = posZ[j];   /* posZ, pos in cartesian coordinates, Angstroms */
    }
  
  NEIGHBORS *NeighborList;
  NeighborList = (NEIGHBORS *) malloc(Natoms*sizeof(NEIGHBORS));
  
  Create_NeighborList((*NList_Cutoff), Natoms, pos, H, NeighborList);
  /*
    Natom is the number of oxygen atoms
    pos is the position of oxygen atoms, cartesian coordinates, size 3*Natoms
    Hmatrix is the box matrix in Angstroms
    Nlist_Cutoff is the cutoff for oxygen sub-lattice in angstroms
  */

  Ql[0] = Steinhardt_Ql(Natoms, 4, ev_pos, ev_force, ev_kappa, cv_pos, *rmin_ql, *rmax_ql, NeighborList, 0, force, pvtens);
  //printf("\n ev = %f, cv = %f, f = %f, en = %f\n", ev_pos[0], cv_pos[0], ev_force[0], Ql[0]);

  Ql[1] = Steinhardt_Ql(Natoms, 6, ev_pos, ev_force, ev_kappa, cv_pos, *rmin_ql, *rmax_ql, NeighborList, 1, force, pvtens);
  //printf("\n ev = %f, cv = %f, f = %f, en = %f\n", ev_pos[1], cv_pos[1], ev_force[1], Ql[1]);

  for (i = 0; i < Natoms; i++)
    {
      j = 9*i;

      // Force is 3*np, while force is 3*Natoms
      Force[j+0] = force[3*i+0]; /* update forces on oxygen atoms */
      Force[j+1] = force[3*i+1];
      Force[j+2] = force[3*i+2];
    }

  /* Stress is 3x3 while pvtens is 9x1 */
  /*
  k = 0;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      Stress[i][j] = pvtens[k++];
  */
  
  // erase table memory
  free(NeighborList);   
  free(pos);
  free(force);
}
