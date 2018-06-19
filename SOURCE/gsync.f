      subroutine gsync()
      use multibead, only: comm_bead
c*********************************************************************
c     
c     barrier / synchronization routine
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith
c
c     wl
c     2001/08/31 11:13:46
c     1.5
c     Exp
c
c*********************************************************************

#include "comms.inc"

#ifndef SERIAL

#ifdef SHMEM

      integer ibar,barrier
      ibar = barrier()

#endif
#ifdef SGISHMEM

       call shmem_barrier_all()

#endif
#if MPI
#ifdef MPIU       
#define MPI_BARRIER MPI_BARRIER_
#endif

      call  MPI_BARRIER(comm_bead,ierr)

#endif

#endif

      return
      end
