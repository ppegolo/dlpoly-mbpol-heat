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
c     systolic pulse version. T3D t.forester sept 1994
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2001/08/31 11:13:48
c     1.4
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
      
#ifdef VAMPIR
      call VTBEGIN(150, ierr)
#endif

#ifndef SERIAL

      safe =.true.
c     
c     load up buffers
      
      j=1
      l=1
      do i=1,natms
        
        if(lstme(l).eq.i)then
          
          l=l+1
          buffer(j+1)=dble(i)
          buffer(j+2)=xxx(i)
          buffer(j+3)=yyy(i)
          buffer(j+4)=zzz(i)
          j=j+4
          
        endif
        
      enddo
c     
c     length of message
      
      buffer(1) = dble(j)
c     
c     array position for incoming messages
      
      mxb = mxbuff/2
c     
c     load initial transfer buffer
      
      call gsync()
c     
c     identity of neighbour node for systolic transfer
      
      jdnode=mod(idnode+1,mxnode)
      
      do k=1,mxnode-1
c     
c     identity of node of origin of incoming data
        
        kdnode=mod(idnode+mxnode-k,mxnode)
c     
c     out going message size
        
        nout = nint(buffer(1))
        
#ifdef INTEL

        msg=irecv(Merge1_tag+k,nin,Ilen)
        call csend(Merge1_tag+k,nout,Ilen,jdnode,0)
        call msgwait(msg)
        
        msg=irecv(Merge1_tag+k,buffer(mxb),Dlen*nin)
        call csend(Merge1_tag+k,buffer(1),Dlen*nout,jdnode,0)
        call msgwait(msg)
        
#endif
#if SHMEM
        
        ibar=barrier()
        call shmem_put(buffer(mxb), buffer(1), nout, jdnode)
        call shmem_udcflush()
        ibar=barrier()
        nin = nint(buffer(mxb))
        
#endif
#if MPI
        call MPI_IRECV(nin,1,MPI_INTEGER,
     x    MPI_ANY_SOURCE,Merge1_tag+k,comm_bead,request,ierr)
        
        call MPI_SEND(nout,1,MPI_INTEGER,jdnode,
     x    Merge1_tag+k,comm_bead,ierr)

        call MPI_WAIT(request,status,ierr)
        
        
        call MPI_IRECV(buffer(mxb),nin,MPI_DOUBLE_PRECISION,
     x    MPI_ANY_SOURCE,Merge1_tag+k,comm_bead,request,ierr)

        call MPI_SEND(buffer(1),nout,MPI_DOUBLE_PRECISION,jdnode,
     x    Merge1_tag+k,comm_bead,ierr)
        
        call MPI_WAIT(request,status,ierr)
        
#endif
c     
c     check buffer array not exceeded
        
        if(nin.gt.mxbuff-mxb) safe =.false.
c     
c     position of first data element in incoming array
        
        nin1 = (nin-1)/4
        j = mxb+1
        
        do j1=1,nin1
          
          i = nint(buffer(j))
          xxx(i)=buffer(j+1)
          yyy(i)=buffer(j+2)
          zzz(i)=buffer(j+3)
          j=j+4
          
        enddo
c     
c     shift new data to start of buffer
        
        do i=1,nin
          
          buffer(i)=buffer(mxb-1+i)
          
        enddo
        
      enddo
c     
c     global check 
      
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,47)
#endif
#ifdef VAMPIR
      call VTEND(150, ierr)
#endif
      return
      end





