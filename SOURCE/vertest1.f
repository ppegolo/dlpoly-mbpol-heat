      subroutine vertest1
     x  (newlst,idnode,mxnode,imcon,natms,delr,cell,
     x  xxx,yyy,zzz,xdf,ydf,zdf,xold,yold,zold)
      
c*********************************************************************
c     
c     DL_POLY subroutime to test for updating of Verlet list
c     replicated data version
c     
c     based on atomic displacements
c     
c     copyright daresbury laboratory 1993
c     author -       t. forester may 1993
c     
c     wl
c     2000/01/18 14:06:00
c     1.4
c     Exp
c     
c*********************************************************************
      
#include "dl_params.inc"
      
      logical newlst,newjob
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension xold(msatms),yold(msatms),zold(msatms)
      dimension cell(9)

      save newjob

      data newjob/.true./

#ifdef VAMPIR
      call VTBEGIN(73, ierr)
#endif
      if((natms+mxnode-1)/mxnode.gt.msatms) call error(idnode,112)
c     
c     set up initial arrays 

      if(newjob) then

        j = 0
        do i = idnode+1,natms,mxnode

          j = j + 1
          xold(j) = xxx(i)
          yold(j) = yyy(i)
          zold(j) = zzz(i)

        enddo

        newjob = .false.
        newlst = .true.

      else
c     
c     integrate velocities 
        j = 0
        do i = idnode+1,natms,mxnode

          j = j+1
          xdf(j) = xold(j)-xxx(i)
          ydf(j) = yold(j)-yyy(i)
          zdf(j) = zold(j)-zzz(i)

        enddo

        call images(imcon,0,1,j,cell,xdf,ydf,zdf)
c     
c     maximum displacement 

        rmax = (delr/2.d0)**2
c     
c     test atomic displacements

        moved = 0

        do k = 1,j
          
          dr = (xdf(k)**2+ydf(k)**2+zdf(k)**2)
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

          j = 0
          do i = idnode+1,natms,mxnode

            j = j + 1
            xold(j) = xxx(i)
            yold(j) = yyy(i)
            zold(j) = zzz(i)

          enddo
          
        endif

      endif

#ifdef VAMPIR
      call VTEND(73, ierr)
#endif
      return
      end
