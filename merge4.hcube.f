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
c     hypercube version - w.smith june 1998
c
c     wl
c     2001/08/31 11:13:49
c     1.5
c     Exp
c
c*********************************************************************

#include "dl_params.inc"
#include "comms.inc"

#ifdef MPI

      integer status(MPI_STATUS_SIZE), request
#ifdef MPIU       
#define MPI_SEND MPI_SEND_
#define MPI_IRECV MPI_IRECV_
#define MPI_WAIT MPI_WAIT_
#endif

#endif
#if SHMEM
      integer barrier
#endif
      real*8 q0(ngrp),q1(ngrp),q2(ngrp),q3(ngrp),buffer(nbuff)

#ifndef SERIAL
#ifdef VAMPIR
      call VTBEGIN(164, ierr)
#endif
c
c     check that buffer is large enough

      nsize=(ngrp+mxnode-1)/mxnode
      if(nbuff.lt.8*nsize)call error(idnode,47)
c
c     load initial transfer buffer

      igrp1 = (idnode*ngrp)/mxnode
      igrp2 = ((idnode+1)*ngrp)/mxnode

#ifdef SHMEM

      j=4*igrp1
      igrp = igrp2-igrp1
      call scopy(igrp,q0(igrp1+1),1,buffer(j+1),4)
      call scopy(igrp,q1(igrp1+1),1,buffer(j+2),4)
      call scopy(igrp,q2(igrp1+1),1,buffer(j+3),4)
      call scopy(igrp,q3(igrp1+1),1,buffer(j+4),4)

#elif SGISHMEM

      j=4*igrp1
      igrp = igrp2-igrp1
      call dcopy(igrp,q0(igrp1+1),1,buffer(j+1),4)
      call dcopy(igrp,q1(igrp1+1),1,buffer(j+2),4)
      call dcopy(igrp,q2(igrp1+1),1,buffer(j+3),4)
      call dcopy(igrp,q3(igrp1+1),1,buffer(j+4),4)

#else

      j=4*igrp1
      do i=igrp1+1,igrp2

         buffer(j+1)=q0(i)
         buffer(j+2)=q1(i)
         buffer(j+3)=q2(i)
         buffer(j+4)=q3(i)
         j=j+4

      enddo

#endif

      call gsync()
c     
c     dimension of hypercube

      ndim=nodedim()
c     
c     pass data to neighbouring nodes
      kk=1

      do k=1,ndim

         nnn=4*(igrp2-igrp1)
         k0=idnode/kk-2*(idnode/(2*kk))
         if(k0.eq.0)then
            
            k1=kk+idnode
            ka=(k1/kk)*kk
            kb=ka+kk
            kgrp1 = (ka*ngrp)/mxnode
            kgrp2 = (kb*ngrp)/mxnode
            mmm=4*(kgrp2-kgrp1)
#if MPI
         call MPI_IRECV(buffer(4*kgrp1+1),mmm,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge_tag+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(4*igrp1+1),nnn,MPI_DOUBLE_PRECISION,k1,
     x        Merge_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif
#ifdef INTEL
            msg1=irecv(Merge_tag1+k,buffer(4*kgrp1+1),Dlen*mmm)
            call csend(Merge_tag2+k,buffer(4*igrp1+1),Dlen*nnn,k1,0)
            call msgwait(msg1)
#endif
#if SHMEM
            ibar=barrier()
            call shmem_put(buffer(4*igrp1+1),buffer(4*igrp1+1),nnn,k1)
            call shmem_udcflush()
            ibar=barrier()
#endif
#if SGISHMEM
            call shmem_barrier_all()
            call shmem_put64(buffer(4*igrp1+1),buffer(4*igrp1+1),nnn,k1)
            call shmem_barrier_all()
#endif

         else                   
            
            k2=idnode-kk
            ka=(k2/kk)*kk
            kb=ka+kk
            kgrp1 = (ka*ngrp)/mxnode
            kgrp2 = (kb*ngrp)/mxnode
            mmm=4*(kgrp2-kgrp1)
#if MPI
         call MPI_IRECV(buffer(4*kgrp1+1),mmm,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge_tag+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(4*igrp1+1),nnn,MPI_DOUBLE_PRECISION,k2,
     x        Merge_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif
#ifdef INTEL
            msg2=irecv(Merge_tag2+k,buffer(4*kgrp1+1),Dlen*mmm)
            call csend(Merge_tag1+k,buffer(4*igrp1+1),Dlen*nnn,k2,0)
            call msgwait(msg2)
#endif
#if SHMEM
            ibar=barrier()
            call shmem_put(buffer(4*igrp1+1),buffer(4*igrp1+1),nnn,k2)
            call shmem_udcflush()
            ibar=barrier()
#endif
#if SGISHMEM
            call shmem_barrier_all()
            call shmem_put64(buffer(4*igrp1+1),buffer(4*igrp1+1),nnn,k2)
            call shmem_barrier_all()
#endif

         endif

         igrp1=min(igrp1,kgrp1)
         igrp2=max(igrp2,kgrp2)
         kk=2*kk
         
      enddo
c
c     unpack buffer array
#ifdef SHMEM

      call scopy(ngrp,buffer(1),4,q0(1),1)
      call scopy(ngrp,buffer(2),4,q1(1),1)
      call scopy(ngrp,buffer(3),4,q2(1),1)
      call scopy(ngrp,buffer(4),4,q3(1),1)

#elif SGISHMEM

      call dcopy(ngrp,buffer(1),4,q0(1),1)
      call dcopy(ngrp,buffer(2),4,q1(1),1)
      call dcopy(ngrp,buffer(3),4,q2(1),1)
      call dcopy(ngrp,buffer(4),4,q3(1),1)

#else

      j=0
      do i=1,ngrp

         q0(i)=buffer(j+1)
         q1(i)=buffer(j+2)
         q2(i)=buffer(j+3)
         q3(i)=buffer(j+4)
         j=j+4

      enddo
#endif
#ifdef VAMPIR
      call VTEND(164, ierr)
#endif
#endif

      return
      end
