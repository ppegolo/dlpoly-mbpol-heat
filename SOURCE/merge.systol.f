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
c     author    - w. smith november 1992.
c     MPI version - t. forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2001/08/31 11:13:48
c     1.5
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

#ifdef VAMPIR
      call VTBEGIN(149, ierr)
#endif

#ifndef SERIAL
c
c     check that buffer is large enough

      nsize=(natms+mxnode-1)/mxnode
      if(nbuff.lt.6*nsize)call error(idnode,47)

c
c     load initial transfer buffer

      j=0

c
c     set up this nodes atoms

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      iatm = iatm2-iatm1+1
	
#ifdef SHMEM

      call scopy(iatm,xxx(iatm1),1,buffer(j+1),3)
      call scopy(iatm,yyy(iatm1),1,buffer(j+2),3)
      call scopy(iatm,zzz(iatm1),1,buffer(j+3),3)

#elif SGISHMEM

      call dcopy(iatm,xxx(iatm1),1,buffer(j+1),3)
      call dcopy(iatm,yyy(iatm1),1,buffer(j+2),3)
      call dcopy(iatm,zzz(iatm1),1,buffer(j+3),3)

#else

      do i=iatm1,iatm2

         buffer(j+1)=xxx(i)
         buffer(j+2)=yyy(i)
         buffer(j+3)=zzz(i)
         j=j+3

      enddo

#endif

      call gsync()
c
c     identity of neighbour node for systolic transfer

      jdnode=mod(idnode+1,mxnode)

      do k=1,mxnode-1

c
c     identity of node of origin of incoming data

         kdnode=mod(idnode+mxnode-k,mxnode)
c
c     identity of incoming  atoms

         katm1 = (kdnode*natms)/mxnode + 1
         katm2 = ((kdnode+1)*natms)/mxnode
         katm = katm2-katm1 + 1
c
c     systolic data pulse to transfer data

#ifdef INTEL

         msg=irecv(Merge_tag+k,buffer(3*nsize+1),3*Dlen*nsize)
         call csend(Merge_tag+k,buffer(1),3*Dlen*nsize,jdnode,0)
         call msgwait(msg)

#endif
#if SHMEM

         ibar=barrier()
         call shmem_put(buffer(3*nsize+1), buffer(1), 3*nsize, jdnode)
         call shmem_udcflush()
         ibar=barrier()

#endif
#if SGISHMEM

         call shmem_put64(buffer(3*nsize+1), buffer(1), 3*nsize, jdnode)
         call shmem_barrier_all()

#endif
#if MPI
         call MPI_IRECV(buffer(3*nsize+1),3*nsize,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge_tag+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(1),3*nsize,MPI_DOUBLE_PRECISION,jdnode,
     x        Merge_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif

c
c     merge the incoming data into current arrays

         j=3*nsize

#ifdef SHMEM

         call scopy(katm,buffer(j+1),3,xxx(katm1),1)
         call scopy(katm,buffer(j+2),3,yyy(katm1),1)
         call scopy(katm,buffer(j+3),3,zzz(katm1),1)

         call scopy(3*nsize,buffer(3*nsize+1),1,buffer(1),1)

#elif SGISHMEM

         call dcopy(katm,buffer(j+1),3,xxx(katm1),1)
         call dcopy(katm,buffer(j+2),3,yyy(katm1),1)
         call dcopy(katm,buffer(j+3),3,zzz(katm1),1)

         call dcopy(3*nsize,buffer(3*nsize+1),1,buffer(1),1)
#else

         do i=katm1,katm2

            xxx(i)=buffer(j+1)
            yyy(i)=buffer(j+2)
            zzz(i)=buffer(j+3)
            j=j+3

         enddo

c
c     shift new data to start of buffer

         do i=1,3*nsize

            buffer(i)=buffer(3*nsize+i)

         enddo

#endif

      enddo

#endif

#ifdef VAMPIR
      call VTEND(149, ierr)
#endif
      return
      end
