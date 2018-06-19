      subroutine centroid_vertest
     x  (newlst,idnode,mxnode,natms,delr,imcon,cell,
     x   xxx,yyy,zzz,xold,yold,zold)
      use multibead
c*********************************************************************
c     
c     DL_POLY subroutine to test for updating of Verlet list
c     replicated data version
c     
c     copyright daresbury laboratory 1993
c     author -       t. forester may 1993
c     
c     wl
c     2000/01/18 14:06:00
c     1.4
c     $Sate: Exp $
c
c     04/12/2006
c     F. Paesani: modified for PIMD/CMD with velocity-verlet
c     
c*********************************************************************
      
#include "dl_params.inc"
#include "mpif.h"
      
      logical newlst,newjob
      
      dimension xxx(natms),yyy(natms),zzz(natms)
      dimension xold(msatms),yold(msatms),zold(msatms)
      dimension cell(9)
      dimension dx(natms),dy(natms),dz(natms)

      save newjob

      data newjob/.true./
#ifdef VAMPIR
      call VTBEGIN(8, ierr)
#endif
      if((natms+mxnode-1)/mxnode.gt.msatms) call error(idnode,112)
c     
c     set up initial arrays 
      iatm1 = (idnode*natms)/mxnode+1
      iatm2 = ((idnode+1)*natms)/mxnode
      if(newjob) then
        j = 0
        do i = iatm1, iatm2
           j = j + 1
           xold(j) = xxx(i)
           yold(j) = yyy(i)
           zold(j) = zzz(i)
        enddo
        newjob = .false.
        newlst = .true.

      else
c     
c        integrate velocities 
         j = 0
         do i = iatm1, iatm2
            j = j+1
            dx(j) = xxx(i)-xold(j)
            dy(j) = yyy(i)-yold(j)
            dz(j) = zzz(i)-zold(j)
        enddo

        call images(imcon,0,1,j,cell,dx,dy,dz)
c     
c       maximum displacement 
        rmax = (delr/2.d0)**2
c     
c       test atomic displacements
        moved = 0
        do k = 1,j
           dr = dx(k)**2+dy(k)**2+dz(k)**2
           if(dr.gt.rmax) moved = moved + 1
        enddo
c     
c       global sum of moved atoms
!        if(mxnode.gt.1) call gisum(moved,1,ibuff)
        call MPI_ALLREDUCE(moved,k,1,MPI_INTEGER,MPI_SUM,
     x                     comm_mb,j)
        moved = k
c     
c       test for new verlet list
        newlst = (moved.ge.2)
c     
c       update stored positions
        if(newlst) then
           j = 0
           do i = iatm1, iatm2
              j = j + 1
              xold(j) = xxx(i)
              yold(j) = yyy(i)
              zold(j) = zzz(i)
           enddo
        endif

      endif
#ifdef VAMPIR
      call VTEND(8, ierr)
#endif
      return
      end
