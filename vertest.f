      subroutine vertest
     x  (newlst,idnode,mxnode,natms,delr,tstep,vxx,vyy,vzz,
     x  xold,yold,zold)
      
c*********************************************************************
c     
c     DL_POLY subroutime to test for updating of Verlet list
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
c*********************************************************************
      
#include "dl_params.inc"
      
      logical newlst,newjob
      
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension xold(msatms),yold(msatms),zold(msatms)

      save newjob

      data newjob/.true./
#ifdef VAMPIR
      call VTBEGIN(8, ierr)
#endif
      if((natms+mxnode-1)/mxnode.gt.msatms) call error(idnode,112)
c     
c     set up initial arrays 

      if(newjob) then

        j = 0
        do i = idnode+1,natms,mxnode

          j = j + 1
          xold(j) = 0.d0
          yold(j) = 0.d0
          zold(j) = 0.d0

        enddo

        newjob = .false.
        newlst = .true.

      else
c     
c     integrate velocities 
        j = 0
        do i = idnode+1,natms,mxnode

          j = j+1
          xold(j) = xold(j)+vxx(i)
          yold(j) = yold(j)+vyy(i)
          zold(j) = zold(j)+vzz(i)

        enddo
c     
c     maximum displacement 

        rmax = (delr/2.d0)**2
c     
c     test atomic displacements

        moved = 0

        do k = 1,j
          
          dr = tstep**2*(xold(k)**2+yold(k)**2+zold(k)**2)
          if(dr.gt.rmax) moved = moved + 1
          
        enddo
c     
c     global sum of moved atoms

        if(mxnode.gt.1) call gisum(moved,1,ibuff)
c     
c     test for new verlet list

        newlst = (moved.ge.2)

c     
c     update stored positions

        if(newlst) then

          do k = 1,j

            xold(k) = 0.d0
            yold(k) = 0.d0
            zold(k) = 0.d0
            
          enddo
          
        endif

      endif
#ifdef VAMPIR
      call VTEND(8, ierr)
#endif
      return
      end
