      subroutine shlmerge
     x     (idnode,mxnode,ntshl,listshl,xxx,yyy,zzz,buffer)
      use multibead, only: comm_bead
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for merging core-shell velocity data
c     to restore data replication on all nodes
c     
c     copyright - daresbury laboratory 1993
c     author    - w. smith february 1993
c     MPI version - w. smith june 1995
c     CPP version - w. smith june 1995
c     
c     wl
c     2001/08/31 11:13:51
c     1.8
c     Exp
c
c***********************************************************************
c     
     
#include "dl_params.inc"
#include "comms.inc"
      
      dimension listshl(mxshl,3)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension buffer(mxbuff)

#ifdef MPI
      integer status(MPI_STATUS_SIZE), request
#ifdef MPIU       
#define MPI_SEND MPI_SEND_
#define MPI_IRECV MPI_IRECV_
#define MPI_WAIT MPI_WAIT_
#endif

#endif
#if  SHMEM

       integer barrier

#endif

#ifndef SERIAL
#ifdef VAMPIR
      call VTBEGIN(124, ierr)
#endif
c     
c     check that buffer is large enough
      
      nsize=8*((ntshl+mxnode-1)/mxnode)
      
      if(mxbuff.lt.2*nsize)call error(idnode,425)

c
c     block indices

      ishl1 = (idnode*ntshl)/mxnode+1
      ishl2 = ((idnode+1)*ntshl)/mxnode

c     
c     load initial transfer buffer
      
      n=0
      m=0
      
      do k=ishl1,ishl2

        m=m+1
        
c     
c     indices of core and shell
        
        i=listshl(m,2)
        j=listshl(m,3)
        buffer(n+1)=dble(i)
        buffer(n+2)=dble(j)
        buffer(n+3)=xxx(i)
        buffer(n+4)=yyy(i)
        buffer(n+5)=zzz(i)
        buffer(n+6)=xxx(j)
        buffer(n+7)=yyy(j)
        buffer(n+8)=zzz(j)
        n=n+8

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
c     systolic data pulse to transfer data
        
#if SHMEM

        ibar=barrier()
        call shmem_put(buffer(nsize+1), buffer(1), nsize, jdnode)
        call shmem_udcflush()
        ibar=barrier()

#endif
#if SGISHMEM

        call shmem_barrier_all()
        call shmem_put64(buffer(nsize+1), buffer(1), nsize, jdnode)
        call shmem_barrier_all()

#endif
#if MPI
        call MPI_IRECV(buffer(nsize+1),nsize,MPI_DOUBLE_PRECISION,
     x    MPI_ANY_SOURCE,Shell_tag+k,comm_bead,request,ierr)

        call MPI_SEND(buffer(1),nsize,MPI_DOUBLE_PRECISION,jdnode,
     x    Shell_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif
#ifdef INTEL

        msg=irecv(Shell_tag+k,buffer(nsize+1),Dlen*nsize)
        call csend(Shell_tag+k,buffer(1),Dlen*nsize,jdnode,0)
        call msgwait(msg)

#endif
c     
c     merge the incoming data into current arrays
        
        n=nsize

c
c     block indices

        kshl1 = (kdnode*ntshl)/mxnode+1
        kshl2 = ((kdnode+1)*ntshl)/mxnode

        do m=kshl1,kshl2
          
          i=nint(buffer(n+1))
          j=nint(buffer(n+2))
          
          xxx(i)=buffer(n+3)
          yyy(i)=buffer(n+4)
          zzz(i)=buffer(n+5)
          xxx(j)=buffer(n+6)
          yyy(j)=buffer(n+7)
          zzz(j)=buffer(n+8)

          n=n+8
          
        enddo
        
c     
c     shift new data to start of buffer
        
        do i=1,nsize
          
          buffer(i)=buffer(nsize+i)
          
        enddo
        
      enddo

#ifdef VAMPIR
      call VTEND(124, ierr)
#endif
#endif
      
      return
      end
