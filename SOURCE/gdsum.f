      subroutine gdsum(aaa,nnn,bbb)
      use multibead, only: comm_bead, bead_rank
c     
c***********************************************************************
c     
c     dl_poly global summation subroutine for MPI - hypercube assumed
c     double precision version
c     
c     copyright - daresbury laboratory 1995
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

#if defined SHMEM || defined SGISHMEM
 #include "dl_params.inc"
#endif     
#include "comms.inc"


#ifdef SHMEM
      integer barrier
#endif

      real*8 aaa(nnn),bbb(nnn)

#ifdef MPI

      integer status(MPI_STATUS_SIZE)

#ifdef MPIU       
#define MPI_allreduce MPI_allreduce_
#endif
      call MPI_allreduce(aaa,bbb,nnn,MPI_DOUBLE_PRECISION,
     x  MPI_SUM,comm_bead,ierror)

        do i = 1,nnn
          aaa(i) = bbb(i)
        enddo

#else
#ifndef SERIAL
c     
c     identify node

      iii=bead_rank

c     
c     pass data to neighbouring nodes
      kk=1

      do k=1,nodedim()

         k0=iii/kk-2*(iii/(2*kk))
         if(k0.eq.0)then
            
            k1=kk+iii
#ifdef INTEL
            msg1=irecv(Dgsum_tag1+k,bbb,Dlen*nnn)
            call csend(Dgsum_tag2+k,aaa,Dlen*nnn,k1,0)
            call msgwait(msg1)
#endif
#if SHMEM
            ibar=barrier()
            call shmem_put(bbb, aaa, nnn, k1)
            call shmem_udcflush()
            ibar=barrier()
#endif
#if SGISHMEM
            call shmem_barrier_all()
            call shmem_put64(bbb, aaa, nnn, k1)
            call shmem_barrier_all()
#endif
c     
c     add data received to local array
         
           do i=1,nnn
            
              aaa(i)=aaa(i)+bbb(i)
            
           enddo

         else                   
            
            k2=iii-kk
#ifdef INTEL
            msg2=irecv(Dgsum_tag2+k,bbb,Dlen*nnn)
            call csend(Dgsum_tag1+k,aaa,Dlen*nnn,k2,0)
            call msgwait(msg2)
#endif
#if SHMEM
            ibar=barrier()
            call shmem_put(bbb, aaa, nnn, k2)
            call shmem_udcflush()
            ibar=barrier()
#endif
#if SGISHMEM
            call shmem_barrier_all()
            call shmem_put64(bbb, aaa, nnn, k2)
            call shmem_barrier_all()
#endif
c     
c     add data received to local array

           do i=1,nnn
            
              aaa(i)=bbb(i)+aaa(i)

           enddo
            
         endif
         kk=2*kk
         
      enddo
#endif
#endif
      return
      end






