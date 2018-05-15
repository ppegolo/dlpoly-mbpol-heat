void Create_NeighborList(double nlist_cutoff, int natoms, 
    double *pos, double H[3][3], NEIGHBORS *NList)
{

  /*-------------------------------------------------------*/
  /* brute force Neighborlist calculation for ab initio MD */
  /*-------------------------------------------------------*/

  /* 
     1. NLIST_CUTOFF: set this to be more than 3rd NN
     2. natoms: no of atoms/ions in the system
     3. pos: pointer to an array of size 3*Natomsx1, positions in Cartesian coords (Angstroms)
     4. MAX_NEIGHBOR_COUNT: is the maximum size of nlist, ndist, etc within NList (about 100)
     5. H: is the H-matrix, box vectors
     6. NList: an array of structures
  */

  int i, j;
  double dx[3], dr, dxc[3], iH[3][3];

  /* constants */

  /* calculate the inverse of H-matrix */
  GetInvH(H, iH);

  /* initialize maximum number of neighbors of all atoms to zero */
  for (i = 0; i < natoms; i++)
    NList[i].nmax = 0;  

  //printf("the box: (%f, %f, %f,      (%f, %f, %f, \n", H[0][0],H[0][1],H[0][2],iH[0][0],iH[0][1],iH[0][2]);
  //printf("         %f, %f, %f,       %f, %f, %f, \n", H[1][0],H[1][1],H[1][2],iH[1][0],iH[1][1],iH[1][2]);
  //printf("         %f, %f, %f)       %f, %f, %f) \n", H[2][0],H[2][1],H[2][2],iH[2][0],iH[2][1],iH[2][2]);
  
  for (i = 0; i < (natoms-1); i++)
    {
      for (j = (i+1); j < natoms; j++)
        {

	  /* get r_ij */
          V3V3subV3(&pos[3*i], &pos[3*j], dxc);	  /* pos, dxc is in Angstroms */
	  
	  /* change to reduced units */
          Cart2Red(dxc, iH, dx);   /* iH in 1/Angstroms, dx in reduced units */

	  /* check PBC, -0.50<dx<0.50 */
          dx[0] = checkPBC_Shifted (dx[0]);
          dx[1] = checkPBC_Shifted (dx[1]);
          dx[2] = checkPBC_Shifted (dx[2]);
          
	  /* change to cart. coords */
          Red2Cart(dx, H, dxc); /* H, dx in Angstroms */

	  dr = NormV3(dxc);     /* dr in Angstroms */
	  
	  if (dr < nlist_cutoff)
	    {
	      NList[i].nlist[NList[i].nmax] = j;
	      NList[i].ndist[3*NList[i].nmax + 0] = dxc[0];
	      NList[i].ndist[3*NList[i].nmax + 1] = dxc[1];
	      NList[i].ndist[3*NList[i].nmax + 2] = dxc[2];
	      NList[i].nmax++;
	      //printf("%d, %d, (%f, %f, %f)\n", i, j, dxc[0],dxc[1],dxc[2]); 
	    }
	}
      
      j = NList[i].nmax - 1;
      //printf("%d, %d\n", i, j);
      //printf("%d, %d, (%f, %f, %f)\n", i, j, NList[i].ndist[3*j], 
      //             NList[i].ndist[3*j + 1], NList[i].ndist[3*j+ 2]);  
      
      if (NList[i].nmax >= max_ncount) 
	{
	  printf("\nSERIOUS PROBLEM : Need to increase the max_ncount\n\n");
	  printf("Nlist: %d \n",NList[i].nmax);
	}
    }  
  return;
}

