      subroutine merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
      use multibead, only: comm_bead
c     
c*********************************************************************
c     
c     dl_poly subroutine for merging together coordinate arrays
c     across a number of processors during rigid body algorithm
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t.forester  november 1993
c     hypercube version. w.smith feb 1998
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2001/08/31 11:13:48
c     1.6
c     Exp
c     
c*********************************************************************
c     
      
#include "dl_params.inc"
#include "comms.inc"
      
      logical safe
      integer lstme(mxatms)
#ifdef SHMEM
      integer barrier
#endif

      real*8 xxx(natms),yyy(natms),zzz(natms),buffer(mxbuff)
      
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
      call VTBEGIN(121, ierr)
#endif
      if(mxnode.eq.1)return
c
c     load transfer buffer

      j=0
      k=1
      safe=.true.
      mblk=mxbuff/mxnode

      do i=1,natms

        if(lstme(k).eq.i)then

          if(j+6.le.mblk)then

            buffer(j+1)=dble(i)
            buffer(j+2)=xxx(i)
            buffer(j+3)=yyy(i)
            buffer(j+4)=zzz(i)
            j=j+4
            k=k+1

          else

            safe=.false.

          endif

        endif

      enddo

      if(safe)buffer(j+1)=dble(0)

      call gstate(safe)
      if(.not.safe)call error(idnode,42)

      mblk=j+1
      call gimax(mblk,1,buffer(j+2))
      iblk=idnode*mblk
      nbuf=mxnode*mblk

      do i=1,j+1

        buffer(iblk+i)=buffer(i)

      enddo
c     
c     pass data to neighbouring nodes
      kk=1
      ip=mblk

      do k=1,nodedim()

        kp=idnode/kk
        k0=idnode/kk-2*((idnode/kk)/2)
        if(k0.eq.0)then
          
          k1=kk+idnode
          ka=kp*ip
          kb=ka+ip
          kc=kb+ip
          ka=min(nbuf,ka)
          kb=min(nbuf,kb)
          kc=min(nbuf,kc)
#if MPI
          call MPI_IRECV(buffer(kb+1),(kc-kb),MPI_DOUBLE_PRECISION,
     x      MPI_ANY_SOURCE,Merge_tag+k,comm_bead,request,ierr)

          call MPI_SEND(buffer(ka+1),ip,MPI_DOUBLE_PRECISION,k1,
     x      Merge_tag+k,comm_bead,ierr)

          call MPI_WAIT(request,status,ierr)
#endif
#ifdef INTEL
          msg1=irecv(Merge_tag1+k,buffer(kb+1),Dlen*(kc-kb))
          call csend(Merge_tag2+k,buffer(ka+1),Dlen*ip,k1,0)
          call msgwait(msg1)
#endif
#if SHMEM
          ibar=barrier()
          call shmem_put(buffer(ka+1),buffer(ka+1),ip,k1)
          call shmem_udcflush()
          ibar=barrier()
#endif
#if SGISHMEM
          call shmem_barrier_all()
          call shmem_put64(buffer(ka+1),buffer(ka+1),ip,k1)
          call shmem_barrier_all()
#endif

        else                   
          
          k2=idnode-kk
          kb=kp*ip
          ka=kb-ip
          kc=kb+ip
          ka=min(nbuf,ka)
          kb=min(nbuf,kb)
          kc=min(nbuf,kc)
#if MPI
          call MPI_IRECV(buffer(ka+1),ip,MPI_DOUBLE_PRECISION,
     x      MPI_ANY_SOURCE,Merge_tag+k,comm_bead,request,ierr)

          call MPI_SEND(buffer(kb+1),(kc-kb),MPI_DOUBLE_PRECISION,k2,
     x      Merge_tag+k,comm_bead,ierr)

          call MPI_WAIT(request,status,ierr)
#endif
#ifdef INTEL
          msg2=irecv(Merge_tag2+k,buffer(ka+1),Dlen*ip)
          call csend(Merge_tag1+k,buffer(kb+1),Dlen*(kc-kb),k2,0)
          call msgwait(msg2)
#endif
#if SHMEM
          ibar=barrier()
          call shmem_put(buffer(kb+1),buffer(kb+1),(kc-kb),k2)
          call shmem_udcflush()
          ibar=barrier()
#endif
#if SGISHMEM
          call shmem_barrier_all()
          call shmem_put64(buffer(kb+1),buffer(kb+1),(kc-kb),k2)
          call shmem_barrier_all()
#endif
        endif

        kk=2*kk
        ip=2*ip
        
      enddo

c     unpack buffer array

      do k=0,mxnode-1

        safe=.true.
        iblk=k*mblk
        jblk=iblk+mblk
        
        do j=iblk,jblk,4
          
          if(safe)then
            
            i=nint(buffer(j+1))
            if(i.gt.0)then
              
              xxx(i)=buffer(j+2)
              yyy(i)=buffer(j+3)
              zzz(i)=buffer(j+4)
              
            else
              
              safe=.false.
              
            endif
            
          endif
          
        enddo

      enddo
#ifdef VAMPIR
      call VTEND(121, ierr)
#endif
#endif

      return
      end
