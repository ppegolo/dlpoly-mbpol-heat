      subroutine csend(msgtag,buf,length,pe,idum)
      use multibead, only: comm_bead
c*********************************************************************
c
c     Intel-like  csend (double precision)
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2001/08/31 11:13:43
c     1.3
c     Exp
c
c*********************************************************************

      implicit none

#include "comms.inc"

      integer msgtag,length,pe,idum
#ifndef SERIAL

      integer len1,ierr
      real*8 buf(*)
#ifdef MPIU       
#define MPI_send MPI_send_
#endif

      len1 = length/Dlen
      call MPI_send(buf,len1,MPI_DOUBLE_PRECISION,pe,msgtag,
     x     comm_bead,ierr)

#else
      integer buf(*)
#endif
      return
      end
