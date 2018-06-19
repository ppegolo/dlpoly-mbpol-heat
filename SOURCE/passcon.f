      subroutine passcon
     x     (lshmov,idnode,mxnode,natms,nscons,lashap,lishap,listme,
     x     listin,listot,listcon,lstfrz)
      use multibead, only: comm_bead
      
c     
c*********************************************************************
c     
c     dl_poly subroutine for passing information about bond 
c     constraints between nodes
c     
c     parallel replicated data version assuming direct node-node
c     connection (i.e. this version may be intel specific)
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith august 1992.
c     MPI version t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2001/08/31 11:13:51
c     1.6
c     Exp
c
c***********************************************************************
c     

#include "dl_params.inc"
#include "comms.inc"

      logical safe,lshmov
      
      dimension listme(mxatms),listin(mxatms),listot(mxatms)
      dimension lishap(mxlshp),lashap(mxproc),listcon(mxcons,3)
      dimension lstfrz(mxatms)

#ifdef MPI
      integer status(MPI_STATUS_SIZE), request
#ifdef MPIU       
#define MPI_SEND MPI_SEND_
#define MPI_IRECV MPI_IRECV_
#define MPI_WAIT MPI_WAIT_
#endif

#endif
#if  defined SHMEM || defined SGISHMEM
      integer status(MPI_STATUS_SIZE)
#endif
      
#ifdef VAMPIR
      call VTBEGIN(153, ierr)
#endif
      if(mxproc.lt.mxnode)call error(idnode,102)
      
      safe=.true.

      do i=1,natms
         
         listme(i)=0
         
      enddo
      
      do k=1,nscons
         
         i=listcon(k,2)
         j=listcon(k,3)
         listme(i)=listme(i)+1
         listme(j)=listme(j)+1
         
      enddo
      
#ifndef SERIAL

      if(mxnode.gt.1)then
         
         j=0
         call gsync()
         do k=1,mxnode-1
            
            jdnode=mod(idnode+mxnode-k,mxnode)
#ifdef INTEL
            msg=irecv(Passcon_tag+k,listin,4*natms)
            call csend(Passcon_tag+k,listme,4*natms,jdnode,0)
            call msgwait(msg)
#endif
#ifdef MPI
            call MPI_IRECV(listin,natms,MPI_INTEGER,
     x        MPI_ANY_SOURCE,Passcon_tag+k,comm_bead,request,ierr)
            
         call MPI_SEND(listme,natms,MPI_INTEGER,jdnode,
     x           Passcon_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif
#if  defined SHMEM || defined SGISHMEM
         call MPI_IRECV(listin,natms,MPI_INTEGER,
     x        MPI_ANY_SOURCE,Passcon_tag+k,comm_bead,request,ierr)

         call MPI_SEND(listme,natms,MPI_INTEGER,jdnode,
     x           Passcon_tag+k,comm_bead,ierr)

         call MPI_WAIT(request,status,ierr)
#endif
            do i=1,natms
               
               if((listme(i).gt.0).and.(listin(i).gt.0.and.
     x              lstfrz(i).eq.0))then
                  
                  j=j+1
                  if(j.gt.mxlshp)then

                     safe=.false.

                  else

                     lishap(j)=i

                  endif
                  
               endif
               
            enddo
            
            lashap(k)=j
            
         enddo
         
      endif

c
c     check for global error condition

      if(mxnode.gt.1) call gstate(safe)

      if(.not.safe)call error(idnode,103)

      if(mxnode.gt.1) then
         call gisum(j,1,idum)
         if(idnode.eq.0) write(nrite,'(/,a,14x,i10)')
     x     ' shared atoms from passcon',j/2
         lshmov = (j.gt.0)
      endif

#endif

c     
c     keep record of all atoms subject to constraints
      
      do i=1,natms
         
         if(listme(i).gt.0)then
            
            listot(i)=1
            
         else
            
            listot(i)=0
            
         endif
         
      enddo
      
      if(mxnode.gt.1)call gisum(listot,natms,listin)

#ifdef VAMPIR
      call VTEND(153, ierr)
#endif
      return
      end
