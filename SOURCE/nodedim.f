      integer function nodedim()
      use multibead, only: comm_bead
c*********************************************************************
c
c     calculate dimension of hypercube
c
c     MPI version - t.forester may 1995
c
c     wl
c     2001/08/31 11:13:50
c     1.5
c     Exp
c
c*********************************************************************

      implicit none
#include "comms.inc"
#ifdef SERIAL
      nodedim = 0
#endif
#ifdef MPI
      integer i,n,ierr,mxnode

#ifdef MPIU       
#define MPI_COMM_SIZE MPI_COMM_SIZE_
#endif
      call MPI_COMM_SIZE(comm_bead, mxnode ,ierr)
      n=1
      nodedim = -1
      do i=0,16

         if(n.eq.mxnode)nodedim=i
         n=2*n

      enddo
c     mpi version does not require hypercube architecture
#endif
#if  defined SHMEM || defined SGISHMEM
      integer i,n,ierr,mxnode

      call MPI_COMM_SIZE(comm_bead, mxnode ,ierr)
      n=1
      nodedim = -1
      do i=0,16

         if(n.eq.mxnode)nodedim=i
         n=2*n

      enddo
      if(nodedim.lt.0)
     x     call error(0,2)

#endif
      return
      end
