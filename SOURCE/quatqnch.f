      subroutine quatqnch
     x  (idnode,imcon,mxnode,natms,ngrp,lstgtp,lstrgd,numgsit,
     x  lstme,gxx,gyy,gzz,buffer,cell,xxt,yyt,zzt,gcmx,gcmy,
     x  gcmz,gmass,gvxx,gvyy,gvzz,q0,q1,q2,q3,omx,omy,omz,
     x  rotinx,rotiny,rotinz,vxx,vyy,vzz,weight,xxx,yyy,zzz)

c***********************************************************************
c     
c     dlpoly routine to convert atomic velocities to rigid body 
c     c.o.m. and angular velocity
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1993.
c     author   - t.forester nov 1993.
c     author   - t.forester dec 1994 : block data.
c     
c     wl
c     2000/01/18 14:05:53
c     1.3
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension gvxx(mxgrp),gvyy(mxgrp),gvzz(mxgrp)
      dimension gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension weight(mxatms),rot(9),buffer(mxbuff)
      dimension gmass(mxungp),omx(mxgrp),omy(mxgrp),omz(mxgrp)
      dimension lstgtp(mxgrp),lstrgd(mxgatm),numgsit(mxungp)
      dimension gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp)
      dimension lstme(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms),cell(9)
#ifdef VAMPIR
      call VTBEGIN(59, ierr)
#endif
c     
c     block indices for groups

      igrp1 = (idnode*ngrp)/mxnode + 1
      igrp2 = ((idnode+1)*ngrp)/mxnode

c     
c     translate atomic velocites to com velocity + angular velocity

      jr=0
      do ig = igrp1,igrp2

        gvxx(ig) =0.d0
        gvyy(ig) =0.d0
        gvzz(ig) =0.d0
        omx(ig) = 0.d0
        omy(ig) = 0.d0
        omz(ig) = 0.d0

        id = lstgtp(ig)

        do j = 1,numgsit(id)

          jr =jr+1
          i =lstrgd(jr)
c     
c     centre of mass momentum

          gvxx(ig) = gvxx(ig) + weight(i)*vxx(i)
          gvyy(ig) = gvyy(ig) + weight(i)*vyy(i)
          gvzz(ig) = gvzz(ig) + weight(i)*vzz(i)
c     
c     distance to c.o.m of molecule

          xxt(jr) = xxx(i) - gcmx(ig)
          yyt(jr) = yyy(i) - gcmy(ig)
          zzt(jr) = zzz(i) - gcmz(ig)

        enddo
c     
c     centre of mass velocity

        gvxx(ig) = gvxx(ig)/gmass(id)
        gvyy(ig) = gvyy(ig)/gmass(id)
        gvzz(ig) = gvzz(ig)/gmass(id)

      enddo

      call images(imcon,0,1,jr,cell,xxt,yyt,zzt)

      jr = 0
      do ig = igrp1,igrp2
c     
c     rotational matrix

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
c     angular momentum accumulators

        wxx = 0.d0
        wyy = 0.d0
        wzz = 0.d0

        id = lstgtp(ig)

        do j = 1,numgsit(id)

          jr =jr+1
          i =lstrgd(jr)

          wxx = wxx + weight(i)*(yyt(jr)*vzz(i) - zzt(jr)*vyy(i))
          wyy = wyy + weight(i)*(zzt(jr)*vxx(i) - xxt(jr)*vzz(i))
          wzz = wzz + weight(i)*(xxt(jr)*vyy(i) - yyt(jr)*vxx(i))

        enddo
c     
c     angular velocity in body fixed frame

        omx(ig) = (rot(1)*wxx+rot(4)*wyy+rot(7)*wzz)*rotinx(id,2)
        omy(ig) = (rot(2)*wxx+rot(5)*wyy+rot(8)*wzz)*rotiny(id,2)
        omz(ig) = (rot(3)*wxx+rot(6)*wyy+rot(9)*wzz)*rotinz(id,2)

        jr = jr -numgsit(id)
        do j = 1,numgsit(id)
          
          jr = jr +1
          i = lstrgd(jr)
          
c     
c     site velocity in body frame 

          wxx = omy(ig)*gzz(id,j) - omz(ig)*gyy(id,j)
          wyy = omz(ig)*gxx(id,j) - omx(ig)*gzz(id,j)
          wzz = omx(ig)*gyy(id,j) - omy(ig)*gxx(id,j)
c     
c     new atomic velocites in lab frame

          vxx(i) = rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+gvxx(ig)
          vyy(i) = rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+gvyy(ig)
          vzz(i) = rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+gvzz(ig)

        enddo

      enddo

      if(mxnode.gt.1) then

        nbuff = mxbuff
        call merge(idnode,mxnode,ngrp,nbuff,gvxx,gvyy,gvzz,buffer)
        call merge(idnode,mxnode,ngrp,nbuff,omx,omy,omz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

      endif

#ifdef VAMPIR
      call VTEND(59, ierr)
#endif
      return
      end
