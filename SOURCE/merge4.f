      subroutine merge4(idnode,mxnode,ngrp,nbuff,q0,q1,q2,q3,buffer)
      use multibead, only: comm_bead
c*********************************************************************
c     
c     dl_poly subroutine for merging coordinate arrays across
c     a number of processors
c     
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1994
c     author    - t.forester  february 1994
c     T3D version - sept 1994 t.forester
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2001/08/31 11:13:49
c     1.8
c     Exp
c
c*********************************************************************

#include "dl_params.inc"
#include "comms.inc"
#ifdef SHMEM
      integer barrier
#endif

      real*8 q0(ngrp),q1(ngrp),q2(ngrp),q3(ngrp),buffer(nbuff)

#ifdef MPI

      integer status(MPI_STATUS_SIZE), request
#ifdef MPIU       
#define MPI_SEND MPI_SEND_
#define MPI_IRECV MPI_IRECV_
#define MPI_WAIT MPI_WAIT_
#endif

#endif

#ifndef SERIAL

#ifdef VAMPIR
      call VTBEGIN(122, ierr)
#endif
c
c     check that buffer is large enough

      nsize=(ngrp+mxnode-1)/mxnode
      if(nbuff.lt.8*nsize)call error(idnode,47)

c
c     load initial transfer buffer

      j=0

      igrp1 = (idnode*ngrp)/mxnode+1
      igrp2 = ((idnode+1)*ngrp)/mxnode

      do i=igrp1,igrp2

         buffer(j+1)=q0(i)
         buffer(j+2)=q1(i)
         buffer(j+3)=q2(i)
         buffer(j+4)=q3(i)
         j=j+4

      enddo

      call gsync()

c
c     identity of neighbour node for systolic transfer

      jdnode=mod(idnode+1,mxnode)
     
      do k=1,mxnode-1

c
c     identity of node of origin of incoming data

         kdnode=mod(idnode+mxnode-k,mxnode)
c
c	identity of incoming groups 

         kgrp1 = (kdnode*ngrp)/mxnode+1
         kgrp2 = ((kdnode+1)*ngrp)/mxnode

#ifdef INTEL
        msg=irecv(Merge4_tag+k,buffer(4*nsize+1),4*Dlen*nsize)
        call csend(Merge4_tag+k,buffer(1),4*Dlen*nsize,jdnode,0)
        call msgwait(msg)
#endif
#if SHMEM

         ibar=barrier()
         call shmem_put(buffer(4*nsize+1), buffer(1), 4*nsize, jdnode)
         call shmem_udcflush()
         ibar=barrier()

#endif
#if SGISHMEM

         call shmem_put64(buffer(4*nsize+1), buffer(1), 4*nsize, jdnode)
         call shmem_barrier_all()

#endif
#if MPI
         call MPI_IRECV(buffer(4*nsize+1),4*nsize,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge4_tag+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(1),4*nsize,MPI_DOUBLE_PRECISION,jdnode,
     x        Merge4_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif

c
c     merge the incoming data into current arrays

         j=4*nsize

         do i=kgrp1,kgrp2

            q0(i)=buffer(j+1)
            q1(i)=buffer(j+2)
            q2(i)=buffer(j+3)
            q3(i)=buffer(j+4)
            j=j+4

         enddo

c
c     shift new data to start of buffer

         do i=1,4*nsize

            buffer(i)=buffer(4*nsize+i)

         enddo

      enddo
#ifdef VAMPIR
      call VTEND(122, ierr)
#endif
#endif
      return
      end
