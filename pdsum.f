      subroutine pdsum3(idnode,mxnode,natms,xxx,yyy,zzz,nbuff,buffer)
      use multibead, only: comm_bead
c
c*********************************************************************
c

#ifdef MPI
#include "mpif.h"
      integer status(MPI_STATUS_SIZE), request
#endif

      real*8 xxx(natms),yyy(natms),zzz(natms),buffer(nbuff)

#ifndef SERIAL

!     check that buffer is large enough
      nsize=(natms+mxnode-1)/mxnode
      if(nbuff.lt.3*(natms+nsize))call error(idnode,47)

      j=0
      do i=1,natms
         buffer(j+1)=xxx(i)
         buffer(j+2)=yyy(i)
         buffer(j+3)=zzz(i)
         j=j+3
      enddo

!     set up this nodes atoms

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      call gsync()

      do k=1,mxnode-1

!        identity of destination node
         kdnode=mod(idnode+mxnode-k,mxnode)

!        outgoing atoms
         katm1 = (kdnode*natms)/mxnode + 1
!         katm2 = ((kdnode+1)*natms)/mxnode

         call MPI_IRECV(buffer(3*natms+1),3*nsize,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,27606+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(3*katm1-2),3*nsize,MPI_DOUBLE_PRECISION,
     x        kdnode, 27606+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
!     merge the incoming data into current arrays

         j=3*natms

         do i=iatm1,iatm2
            xxx(i)=xxx(i)+buffer(j+1)
            yyy(i)=yyy(i)+buffer(j+2)
            zzz(i)=zzz(i)+buffer(j+3)
            j=j+3
         enddo
      enddo

#endif
      return
      end

      subroutine pdsum6(idnode,mxnode,natms,xxx,yyy,zzz,
     x                  uuu,vvv,ttt,nbuff,buffer)
      use multibead, only: comm_bead
c
c*********************************************************************
c

#ifdef MPI
#include "mpif.h"
      integer status(MPI_STATUS_SIZE), request
#endif

      real*8 xxx(natms),yyy(natms),zzz(natms),
     x       uuu(natms),vvv(natms),ttt(natms), buffer(nbuff)

#ifndef SERIAL

!     check that buffer is large enough
      nsize=(natms+mxnode-1)/mxnode
      if(nbuff.lt.6*(natms+nsize))call error(idnode,47)

      j=0
      do i=1,natms
         buffer(j+1)=xxx(i)
         buffer(j+2)=yyy(i)
         buffer(j+3)=zzz(i)
         buffer(j+4)=uuu(i)
         buffer(j+5)=vvv(i)
         buffer(j+6)=ttt(i)
         j=j+6
      enddo

!     set up this nodes atoms

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      call gsync()

      do k=1,mxnode-1

!        identity of destination node
         kdnode=mod(idnode+mxnode-k,mxnode)

!        outgoing atoms
         katm1 = (kdnode*natms)/mxnode + 1
!         katm2 = ((kdnode+1)*natms)/mxnode

         call MPI_IRECV(buffer(6*natms+1),6*nsize,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,27606+k,comm_bead,request,ierr)

         call MPI_SEND(buffer(6*katm1-5),6*nsize,MPI_DOUBLE_PRECISION,
     x        kdnode, 27606+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
!     merge the incoming data into current arrays

         j=6*natms

         do i=iatm1,iatm2
            xxx(i)=xxx(i)+buffer(j+1)
            yyy(i)=yyy(i)+buffer(j+2)
            zzz(i)=zzz(i)+buffer(j+3)
            uuu(i)=uuu(i)+buffer(j+4)
            vvv(i)=vvv(i)+buffer(j+5)
            ttt(i)=ttt(i)+buffer(j+6)
            j=j+6
         enddo
      enddo

#endif
      return
      end
