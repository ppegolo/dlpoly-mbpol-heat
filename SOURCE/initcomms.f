      subroutine initcomms()
      
c*********************************************************************
c     
c     communication harness initialisation
c     
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2001/08/31 11:13:47
c     1.6
c     Exp
c     
c*********************************************************************
      
      implicit none
      
#include "comms.inc"
      
#ifndef SERIAL
      
#if defined MPI || defined SHMEM || defined SGISHMEM
      integer ierr

c     required for MPI only
#ifdef MPIU       
#define MPI_init MPI_init_
#endif
      call MPI_init(ierr)
      
#endif
#endif

      return
      end



c PLUMED
      integer*8 function get_comms()

c*********************************************************************
c     
c     dl_poly subroutine for obtaining the communicator
c     this is used by plumed
c     
c     author - G. Tribello 
c
c*********************************************************************

      get_comms=MPI_COMM_WORLD
      end

c PLUMED




