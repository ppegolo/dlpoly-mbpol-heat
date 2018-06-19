      subroutine crecv(msgtag,buf,length)
      use multibead, only: comm_bead
c*********************************************************************
c
c     Intel-like  crecv (double precision)
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

      integer msgtag,length
#ifndef SERIAL

      integer len1,ierr
      integer status(MPI_STATUS_SIZE)
      real*8 buf(*)

      len1 = length/Dlen

#ifdef MPIU       
#define MPI_RECV MPI_RECV_
#endif
      call MPI_RECV(buf,len1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
     x     msgtag,comm_bead,status,ierr)

#else
      integer buf(*)
#endif
      return 
      end
