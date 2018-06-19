      subroutine gimax(aaa,nnn,bbb)
      use multibead, only: comm_bead
c     
c***********************************************************************
c     
c     dl_poly global maximum subroutine for hypercube - MPI version
c     integer version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2001/08/31 11:13:45
c     1.6
c     Exp
c
c***********************************************************************
c     
      
#include "dl_params.inc"
#include "comms.inc"

#ifndef SERIAL

#if defined MPI || defined SHMEM || defined SGISHMEM
      integer status(MPI_STATUS_SIZE)
#ifdef MPIU       
#define MPI_allreduce MPI_allreduce_
#endif
      integer aaa(nnn),bbb(nnn)

      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x  MPI_MAX,comm_bead,ierror)

        do i = 1,nnn
          aaa(i) = bbb(i)
        enddo

#endif

#ifdef INTEL
      integer aaa(nnn),bbb(nnn)
c     
c     identify node
      iii=idnode

c     
c     pass data to neighbouring nodes
      kk=1

      call gsync()

      do k=1,nodedim()

         k0=iii/kk-2*(iii/(2*kk))
         if(k0.eq.0)then
            
            k1=kk+iii
            msg1=irecv(Igsum_tag1+k,bbb,Ilen*nnn)
            call csend(Igsum_tag2+k,aaa,Ilen*nnn,k1,0)
            call msgwait(msg1)

         else                   
            
            k2=iii-kk
            msg2=irecv(Igsum_tag2+k,bbb,Ilen*nnn)
            call csend(Igsum_tag1+k,aaa,Ilen*nnn,k2,0)
            call msgwait(msg2)
            
         endif
c     
c     take maximum
         
         do i=1,nnn
           
           aaa(i)=max(bbb(i),aaa(i))
           
         enddo
         kk=2*kk
         
      enddo
#endif

#endif
      return
      end
