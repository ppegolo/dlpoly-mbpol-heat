      subroutine strucopt
     x  (loptim,idnode,mxnode,imcon,natms,ngrp,ntfree,
     x  lstfre,lstgtp,lstrgd,numgsit,cell,gvxx,gvyy,gvzz,
     x  gcmx,gcmy,gcmz,omx,omy,omz,vxx,vyy,vzz,fxx,fyy,
     x  fzz,xxx,yyy,zzz,xxt,yyt,zzt,q0,q1,q2,q3,rotinx,
     x  rotiny,rotinz)
      
c***********************************************************************
c     
c     dl_poly routine for zero Kelvin temperature optimization
c     if veloc.Force < 0 then velocity is set to zero in 
c     preparation for integration of equations of motion
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1994
c     author t.forester     march 1994
c     amended t.forester    dec 1994 : block data
c     
c     wl
c     2000/01/18 14:05:57
c     1.4
c     Exp
c     
c***********************************************************************
      
#include "dl_params.inc"
      
      logical loptim
      dimension lstfre(mxatms),lstrgd(mxgatm),numgsit(mxungp)
      dimension lstgtp(mxgrp)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension gvxx(mxgrp),gvyy(mxgrp),gvzz(mxgrp)
      dimension gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp)
      dimension omx(mxgrp),omy(mxgrp),omz(mxgrp)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension cell(9),rot(9)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)

#ifdef VAMPIR
      call VTBEGIN(30, ierr)
#endif
      if(loptim) then 
        
        if(ngrp.eq.0) then
c     
c     block indices

          iatm1 = (idnode*natms)/mxnode + 1
          iatm2 = ((idnode+1)*natms)/mxnode
c     
c     take component of veloc in direction of force
          
          do i = iatm1,iatm2

            dot = vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
            if(dot.lt.0.d0) then
              vxx(i) = 0.d0
              vyy(i) = 0.d0
              vzz(i) = 0.d0
            else
c
c     take component of velocity in direction of force

              fsq = (fxx(i)**2+fyy(i)**2+fzz(i)**2)
              fsq = dot/max(1.d-10,fsq)
              vxx(i) = fxx(i)*fsq
              vyy(i) = fyy(i)*fsq
              vzz(i) = fzz(i)*fsq
            endif
            
          enddo
          
        else
c     
c     block indices for groups and free atoms

          igrp1 = (idnode*ngrp)/mxnode + 1
          igrp2 = ((idnode+1)*ngrp)/mxnode

          ifre1 = (idnode*ntfree)/mxnode + 1
          ifre2 = ((idnode+1)*ntfree)/mxnode

          do j = ifre1,ifre2
c     
c     reset atomic velocities 
            
            i = lstfre(j)
            
            dot = vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
            if(dot.lt.0.d0) then
              vxx(i) = 0.d0
              vyy(i) = 0.d0
              vzz(i) = 0.d0
            else
c
c     take component of velocity in direction of force

              fsq = (fxx(i)**2+fyy(i)**2+fzz(i)**2)
              fsq = dot/max(1.d-10,fsq)
              vxx(i) = fxx(i)*fsq
              vyy(i) = fyy(i)*fsq
              vzz(i) = fzz(i)*fsq
            endif
            
          enddo

          jr = 0
          do ig = igrp1,igrp2
c     
c     reset rigid body velocites (linear and angular)
            
            fcomx = 0.d0
            fcomy = 0.d0
            fcomz = 0.d0

            id =lstgtp(ig)
            do j = 1,numgsit(id)
              
              jr = jr+1
              i = lstrgd(jr)
c     
c     forces on com
              
              fcomx = fcomx + fxx(i)
              fcomy = fcomy + fyy(i)
              fcomz = fcomz + fzz(i)
              
            enddo
            
            dot = gvxx(ig)*fcomx+gvyy(ig)*fcomy+gvzz(ig)*fcomz
            if(dot.lt.0.d0) then
              gvxx(ig) = 0.d0
              gvyy(ig) = 0.d0
              gvzz(ig) = 0.d0
            else
c
c     take component of velocity in direction of force

              fsq = (fcomx**2+fcomy**2+fcomz**2)
              fsq = dot/max(1.d-10,fsq)
              gvxx(ig) = fcomx*fsq
              gvyy(ig) = fcomy*fsq
              gvzz(ig) = fcomz*fsq

            endif

          enddo
c     
c     calculate torques
c     
c     site to com distances
          
          jr=0
          do ig = igrp1,igrp2
            
            do j = 1,numgsit(lstgtp(ig))
              
              jr = jr +1
              i = lstrgd(jr)
              
              xxt(jr) = xxx(i) - gcmx(ig)
              yyt(jr) = yyy(i) - gcmy(ig)
              zzt(jr) = zzz(i) - gcmz(ig)
              
            enddo
            
          enddo
c     
c     minimum images
          
          call images(imcon,0,1,jr,cell,xxt,yyt,zzt)
c     
c     torques in lab frame
          
          jr=0
          do ig = igrp1,igrp2
            
            trx = 0.d0
            try = 0.d0
            trz = 0.d0
            
            id = lstgtp(ig)
            do j = 1,numgsit(id)
              
              jr = jr +1
              i = lstrgd(jr)
              
              trx = trx + yyt(jr)*fzz(i) - zzt(jr)*fyy(i)
              try = try + zzt(jr)*fxx(i) - xxt(jr)*fzz(i)
              trz = trz + xxt(jr)*fyy(i) - yyt(jr)*fxx(i)
              
            enddo
            
            rot(1) = q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
            rot(2) = 2.d0*(q1(ig)*q2(ig) - q0(ig)*q3(ig))
            rot(3) = 2.d0*(q1(ig)*q3(ig) + q0(ig)*q2(ig))
            rot(4) = 2.d0*(q1(ig)*q2(ig) + q0(ig)*q3(ig))
            rot(5) = q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
            rot(6) = 2.d0*(q2(ig)*q3(ig) - q0(ig)*q1(ig))
            rot(7) = 2.d0*(q1(ig)*q3(ig) - q0(ig)*q2(ig))
            rot(8) = 2.d0*(q2(ig)*q3(ig) + q0(ig)*q1(ig))
            rot(9) = q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
c     
c     transform to body fixed frame
            
            tax=(trx*rot(1)+try*rot(4)+trz*rot(7))*rotinx(id,2)
            tay=(trx*rot(2)+try*rot(5)+trz*rot(8))*rotiny(id,2)
            taz=(trx*rot(3)+try*rot(6)+trz*rot(9))*rotinz(id,2)
            
            dot = omx(ig)*tax + omy(ig)*tay + omz(ig)*taz
            if(dot.le.0.d0) then
              omx(ig) = 0.d0
              omy(ig) = 0.d0
              omz(ig) = 0.d0
            else
c     
c     take component of velocity in direction of torque
              
              fsq = (tax**2+tay**2+taz**2)
              fsq = dot/max(1.d-10,fsq)
              omx(ig) = tax*fsq
              omy(ig) = tay*fsq
              omz(ig) = taz*fsq
              
            endif
            
          enddo
          
        endif
        
      endif
      
#ifdef VAMPIR
      call VTEND(30, ierr)
#endif
      return
      end

