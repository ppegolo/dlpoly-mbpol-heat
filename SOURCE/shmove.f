      subroutine shmove
     x     (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x      txx,tyy,tzz,buffer)
      use multibead, only: comm_bead
c     
c***********************************************************************
c     
c     dl_poly subroutine for passing coordinate updates between
c     nodes during the shake iteration cycle
c     
c     parallel replicated data algorithm 
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith august 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2001/08/31 11:13:52
c     1.10
c     Exp
c
c***********************************************************************
c     

#include "dl_params.inc"
#include "comms.inc"

      integer idnode, mxnode, natms
      integer lishap(mxlshp),lashap(mxproc)

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

      real*8 xxt(mxatms),yyt(mxatms),zzt(mxatms)
      real*8 txx(mxatms),tyy(mxatms),tzz(mxatms)
      real*8 buffer(mxbuff)

#ifndef SERIAL
#ifdef VAMPIR
      call VTBEGIN(125, ierr)
#endif

c
c     store increments to be transferred

      do i=1,natms

         txx(i)=xxt(i)
         tyy(i)=yyt(i)
         tzz(i)=zzt(i)

      enddo

c
c     transfer coordinate data to all nodes

      call gsync()

      do k=1,mxnode-1

         i=0
         j0=0
         if(k.gt.1)j0=lashap(k-1)

         do j=j0+1,lashap(k)

            buffer(i+1)=txx(lishap(j))
            buffer(i+2)=tyy(lishap(j))
            buffer(i+3)=tzz(lishap(j))
            i=i+3                             
                                               
         enddo

c
c     inter node communication

         k0=0

         if(k+1.lt.mxnode)k0=lashap(mxnode-k-1)
         n=3*(lashap(mxnode-k)-k0)
         jdnode=mod(idnode+k,mxnode)

c
c     check for zero length messages
#ifdef INTEL
         if(n.gt.0) msg=irecv(Shmove_tag+k,buffer(i+1),Dlen*n)
         if(i.gt.0) call csend(Shmove_tag+k,buffer(1),Dlen*i,jdnode,0)
         if(n.gt.0) call msgwait(msg)
#endif
#if SHMEM
         lim=3*mxatms
         ibar=barrier()
         if(i.gt.0) call shmem_put(buffer(lim+1),buffer(1),i,jdnode)
         call shmem_udcflush()
         ibar=barrier()
         if(n.gt.0) call scopy(n,buffer(lim+1),1,buffer(1),1)
         i=0
#endif
#if SGISHMEM
         lim=3*mxatms
         call shmem_barrier_all()
         if(i.gt.0) call shmem_put64(buffer(lim+1),buffer(1),i,jdnode)
         call shmem_barrier_all()
         if(n.gt.0) call dcopy(n,buffer(lim+1),1,buffer(1),1)
         i=0
#endif
#if MPI
         if(n.gt.0) call MPI_IRECV(buffer(i+1),n,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Shmove_tag+k,comm_bead,request,ierr)

         if(i.gt.0) call MPI_SEND(buffer(1),i,MPI_DOUBLE_PRECISION,
     x        jdnode,Shmove_tag+k,comm_bead,ierr)

         if(n.gt.0) call MPI_WAIT(request,status,ierr)
#endif
c
c     consolidate transferred data

         do j=k0+1,lashap(mxnode-k)

            xxt(lishap(j))=xxt(lishap(j))+buffer(i+1)
            yyt(lishap(j))=yyt(lishap(j))+buffer(i+2)
            zzt(lishap(j))=zzt(lishap(j))+buffer(i+3)
            i=i+3

         enddo
         
      enddo

#ifdef VAMPIR
      call VTEND(125, ierr)
#endif
#endif

      return
      end
