      subroutine merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
      use multibead, only: comm_bead
c
c*********************************************************************
c     
c     dl_poly subroutine for merging coordinate arrays across
c     a number of processors
c     
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith november 1992
c     hypercube version - w.smith feb 1998
c     MPI version - t. forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2001/08/31 11:13:47
c     1.6
c     Exp
c
c*********************************************************************
c
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
      real*8 xxx(natms),yyy(natms),zzz(natms),buffer(nbuff)

#ifndef SERIAL
#ifdef VAMPIR
      call VTBEGIN(164, ierr)
#endif
c
c     check that buffer is large enough

      if(nbuff.lt.6*natms)call error(idnode,47)
c
c     set up this nodes atoms

      iatm1 = (idnode*natms)/mxnode
      iatm2 = ((idnode+1)*natms)/mxnode
c
c     load initial transfer buffer
	
#ifdef SHMEM

      j=3*iatm1
      iatm = iatm2-iatm1
      call scopy(iatm,xxx(iatm1+1),1,buffer(j+1),3)
      call scopy(iatm,yyy(iatm1+1),1,buffer(j+2),3)
      call scopy(iatm,zzz(iatm1+1),1,buffer(j+3),3)

#elif SGISHMEM

      j=3*iatm1
      iatm = iatm2-iatm1
      call dcopy(iatm,xxx(iatm1+1),1,buffer(j+1),3)
      call dcopy(iatm,yyy(iatm1+1),1,buffer(j+2),3)
      call dcopy(iatm,zzz(iatm1+1),1,buffer(j+3),3)

#else

      j=3*iatm1
      do i=iatm1+1,iatm2

         buffer(j+1)=xxx(i)
         buffer(j+2)=yyy(i)
         buffer(j+3)=zzz(i)
         j=j+3

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

         nnn=3*(iatm2-iatm1)
         k0=idnode/kk-2*(idnode/(2*kk))
         if(k0.eq.0)then
            
            k1=kk+idnode
            ka=(k1/kk)*kk
            kb=ka+kk
            katm1 = (ka*natms)/mxnode
            katm2 = (kb*natms)/mxnode
            mmm=3*(katm2-katm1)
#if MPI
         call MPI_IRECV(buffer(3*katm1+1),mmm,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge_tag+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(3*iatm1+1),nnn,MPI_DOUBLE_PRECISION,k1,
     x        Merge_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif
#ifdef INTEL
            msg1=irecv(Merge_tag1+k,buffer(3*katm1+1),Dlen*mmm)
            call csend(Merge_tag2+k,buffer(3*iatm1+1),Dlen*nnn,k1,0)
            call msgwait(msg1)
#endif
#if SHMEM
            ibar=barrier()
            call shmem_put(buffer(3*iatm1+1),buffer(3*iatm1+1),nnn,k1)
            call shmem_udcflush()
            ibar=barrier()
#endif
#if SGISHMEM
            call shmem_barrier_all()
            call shmem_put64(buffer(3*iatm1+1),buffer(3*iatm1+1),nnn,k1)
            call shmem_barrier_all()
#endif

         else                   
            
            k2=idnode-kk
            ka=(k2/kk)*kk
            kb=ka+kk
            katm1 = (ka*natms)/mxnode
            katm2 = (kb*natms)/mxnode
            mmm=3*(katm2-katm1)
#if MPI
         call MPI_IRECV(buffer(3*katm1+1),mmm,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge_tag+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(3*iatm1+1),nnn,MPI_DOUBLE_PRECISION,k2,
     x        Merge_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif
#ifdef INTEL
            msg2=irecv(Merge_tag2+k,buffer(3*katm1+1),Dlen*mmm)
            call csend(Merge_tag1+k,buffer(3*iatm1+1),Dlen*nnn,k2,0)
            call msgwait(msg2)
#endif
#if SHMEM
            ibar=barrier()
            call shmem_put(buffer(3*iatm1+1),buffer(3*iatm1+1),nnn,k2)
            call shmem_udcflush()
            ibar=barrier()
#endif
#if SGISHMEM
            call shmem_barrier_all()
            call shmem_put64(buffer(3*iatm1+1),buffer(3*iatm1+1),nnn,k2)
            call shmem_barrier_all()
#endif

         endif

         iatm1=min(iatm1,katm1)
         iatm2=max(iatm2,katm2)
         kk=2*kk
         
      enddo
c
c     unpack buffer array
#ifdef SHMEM

      call scopy(natms,buffer(1),3,xxx(1),1)
      call scopy(natms,buffer(2),3,yyy(1),1)
      call scopy(natms,buffer(3),3,zzz(1),1)

#elif SGISHMEM

      call dcopy(natms,buffer(1),3,xxx(1),1)
      call dcopy(natms,buffer(2),3,yyy(1),1)
      call dcopy(natms,buffer(3),3,zzz(1),1)

#else

      j=0
      do i=1,natms

         xxx(i)=buffer(j+1)
         yyy(i)=buffer(j+2)
         zzz(i)=buffer(j+3)
         j=j+3

      enddo
#endif
#ifdef VAMPIR
      call VTEND(164, ierr)
#endif
#endif

      return
      end
