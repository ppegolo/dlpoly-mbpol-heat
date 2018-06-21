       subroutine qshake
     x  (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x  nscons,tolnce,tstep,vircon,lashap,lishap,listcon,
     x  listme,listot,lstfrz,lstbod,lstgtp,lstcsit,
     x  buffer,cell,dxt,dxx,dyt,dyy,dzt,dzz,prmcon,txx,
     x  tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,fxx,fyy,
     x  fzz,stress,q0,q1,q2,q3,gxx,gyy,gzz,gmass,rotinx,
     x  rotiny,rotinz,redmass,esig1)
c     
c***********************************************************************
c     
c     dl_poly subroutine for appling bond constraint corrections after
c     atomic integration. Assumes rigid bodies connected by constraints
c     If this is not so use rdshake_1 instead
c     Must be used in conjunction with integration algorithms
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester june 1995
c     
c     wl
c     2000/01/18 14:05:53
c     1.5
c     $Sate: Exp $
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe,lshmov,newstep
      
      dimension listme(mxatms)
      dimension lishap(mxlshp)
      dimension lashap(mxproc),listot(mxatms)
      dimension listcon(mxcons,3),buffer(mxbuff)
      dimension lstfrz(mxatms)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension weight(mxatms),prmcon(mxtcon),cell(9)
      dimension stress(9)
      dimension rot(9)
      dimension q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp)
      dimension gxx(mxungp,mxngp),gyy(mxungp,mxngp),gzz(mxungp,mxngp)
      dimension lstbod(mxatms),lstcsit(2*mxcons)
      dimension lstgtp(mxgrp),gmass(mxungp)
      dimension rotinx(mxungp,2),rotiny(mxungp,2),rotinz(mxungp,2)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension redmass(mxcons),esig1(mxcons)
#ifdef VAMPIR
      call VTBEGIN(128, ierr)
#endif
c     
c     constraint virial
      
      vircon=0.d0
#ifdef STRESS
c     
c     accumulators for stress tensor
      strs1 = 0.d0
      strs2 = 0.d0
      strs3 = 0.d0
      strs5 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0
#endif
c     
c     timestep squared
      
      tstep2 = tstep*tstep
c     
c     assume bond vectors dxx,dyy,dzz are input
c     dxx = xxx(i) - xxx(j) etc
      
      safe=.false.
c     
c     one iteration of constraint (shake) algorithm
c     do icyc = 1,mxshak
      
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
c     
c     calculate temporary bond vector
        
        dxt(k)=xxx(i)-xxx(j)
        dyt(k)=yyy(i)-yyy(j)
        dzt(k)=zzz(i)-zzz(j)
        
      enddo
c     
c     periodic boundary condition
      
      call images(imcon,0,1,nscons,cell,dxt,dyt,dzt)
c     
c     calculate maximum error in bondlength
      
      esig=0.d0
      do k=1,nscons
c     
c     set bond parameter
        
        dis  = prmcon(listcon(k,1))
        dis2 = dis*dis
        esig1(k)=0.5d0*(dis2- (dxt(k)**2+dyt(k)**2+dzt(k)**2))/dis2
        esig = max(esig,abs(esig1(k)))
        
      enddo

c     
c     global verification of convergence

      safe=(esig.lt.tolnce)
      
      if(mxnode.gt.1)then
        
        call gstate(safe)
        
      endif
c     
c     terminate iteration if all tolerances satisfied 
      
      if (safe) go to 100 
c     
c     initialise force increment arrays
      
      do i=1,natms
        
        xxt(i)=0.d0
        yyt(i)=0.d0
        zzt(i)=0.d0
        
      enddo
c     
c     calculate constraint forces
      
      ik = 0
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
c
c     assign effective reduced mass

        if(newstep) then
          ig = lstbod(i)

          if(ig.eq.0) then
            amti= 1.d0/weight(i)
          else
          
            ik = ik+1
            id = lstgtp(ig)

            rot(1) = q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
            rot(2) = 2.d0*(q1(ig)*q2(ig) - q0(ig)*q3(ig))
            rot(3) = 2.d0*(q1(ig)*q3(ig) + q0(ig)*q2(ig))
            rot(4) = 2.d0*(q1(ig)*q2(ig) + q0(ig)*q3(ig))
            rot(5) = q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
            rot(6) = 2.d0*(q2(ig)*q3(ig) - q0(ig)*q1(ig))
            rot(7) = 2.d0*(q1(ig)*q3(ig) - q0(ig)*q2(ig))
            rot(8) = 2.d0*(q2(ig)*q3(ig) + q0(ig)*q1(ig))
            rot(9) = q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
          
            jj = lstcsit(ik)

c     
c     site to com in lab frame
            xxa=(gxx(id,jj)*rot(1)+gyy(id,jj)*rot(2)+gzz(id,jj)*rot(3))
            yya=(gxx(id,jj)*rot(4)+gyy(id,jj)*rot(5)+gzz(id,jj)*rot(6))
            zza=(gxx(id,jj)*rot(7)+gyy(id,jj)*rot(8)+gzz(id,jj)*rot(9))
c     
c     find cross product between interatomic vector and vector to com
          
            tax  = yya*dzz(k) - zza*dyy(k)
            tay  = zza*dxx(k) - xxa*dzz(k)
            taz  = xxa*dyy(k) - yya*dxx(k)
c     
c     transform to body fixed frame
          
            trx=(tax*rot(1)+tay*rot(4)+taz*rot(7))*rotinx(id,2)
            try=(tax*rot(2)+tay*rot(5)+taz*rot(8))*rotiny(id,2)
            trz=(tax*rot(3)+tay*rot(6)+taz*rot(9))*rotinz(id,2)
c     
c     direction of induced velocites in body frame
          
            vix = try*gzz(id,jj) - trz*gyy(id,jj)
            viy = trz*gxx(id,jj) - trx*gzz(id,jj)
            viz = trx*gyy(id,jj) - try*gxx(id,jj)
c
c     transform to lab frame

            vxi = vix*rot(1) + viy*rot(2) + viz*rot(3)
            vyi = vix*rot(4) + viy*rot(5) + viz*rot(6)
            vzi = vix*rot(7) + viy*rot(8) + viz*rot(9)
c     
c     find dot product between induced translational and rotational velocities
          
            doti = abs(vxi*dxx(k)+vyi*dyy(k)+vzi*dzz(k))
            doti = doti/dis2
          
            amti = (1.d0/gmass(id) + doti)
            
          endif
        
          ig = lstbod(j)
          if(ig.eq.0) then
          
            amtj= 1.d0/weight(j)

          else
          
            ik = ik+1
            id = lstgtp(ig)
          
            rot(1) = q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
            rot(2) = 2.d0*(q1(ig)*q2(ig) - q0(ig)*q3(ig))
            rot(3) = 2.d0*(q1(ig)*q3(ig) + q0(ig)*q2(ig))
            rot(4) = 2.d0*(q1(ig)*q2(ig) + q0(ig)*q3(ig))
            rot(5) = q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
            rot(6) = 2.d0*(q2(ig)*q3(ig) - q0(ig)*q1(ig))
            rot(7) = 2.d0*(q1(ig)*q3(ig) - q0(ig)*q2(ig))
            rot(8) = 2.d0*(q2(ig)*q3(ig) + q0(ig)*q1(ig))
            rot(9) = q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
          
            jj = lstcsit(ik)
c     
c     site to com in lab frame
            xxa=(gxx(id,jj)*rot(1)+gyy(id,jj)*rot(2)+gzz(id,jj)*rot(3))
            yya=(gxx(id,jj)*rot(4)+gyy(id,jj)*rot(5)+gzz(id,jj)*rot(6))
            zza=(gxx(id,jj)*rot(7)+gyy(id,jj)*rot(8)+gzz(id,jj)*rot(9))
c     
c     find cross product between interatomic vector and vector to com
          
            tax  = yya*dzz(k) - zza*dyy(k)
            tay  = zza*dxx(k) - xxa*dzz(k)
            taz  = xxa*dyy(k) - yya*dxx(k)
c     
c     transform to body fixed frame
          
            trx=(tax*rot(1)+tay*rot(4)+taz*rot(7))*rotinx(id,2)
            try=(tax*rot(2)+tay*rot(5)+taz*rot(8))*rotiny(id,2)
            trz=(tax*rot(3)+tay*rot(6)+taz*rot(9))*rotinz(id,2)
c     
c     direction of induced velocites in body frame
          
            vjx = try*gzz(id,jj) - trz*gyy(id,jj)
            vjy = trz*gxx(id,jj) - trx*gzz(id,jj)
            vjz = trx*gyy(id,jj) - try*gxx(id,jj)
c
c     transform to lab frame

            vxj = vjx*rot(1) + vjy*rot(2) + vjz*rot(3)
            vyj = vjx*rot(4) + vjy*rot(5) + vjz*rot(6)
            vzj = vjx*rot(7) + vjy*rot(8) + vjz*rot(9)
c     
c     find dot product between induced translational and rotational velocities
          
            doti = abs(vxj*dxx(k)+vyj*dyy(k)+vzj*dzz(k))
            doti = doti/dis2

            amtj = (1.d0/gmass(id) + doti)

          endif
        
          if(lstfrz(i).ne.0) then
            amti = 0.d0
          endif

          if(lstfrz(j).ne.0) then 
            amtj = 0.d0
          endif

          redmass(k) = 1.d0/(amti + amtj)/tstep2

        endif
c     
c     constraint force parameter 
        
        gamma = esig1(k)*redmass(k)
c     
c     accumulate bond virial
        
        vircon=vircon-gamma*(dxx(k)**2+dyy(k)**2+dzz(k)**2)
#ifdef STRESS
        strs1 = strs1 + gamma*dxx(k)*dxx(k)
        strs2 = strs2 + gamma*dxx(k)*dyy(k)
        strs3 = strs3 + gamma*dxx(k)*dzz(k)
        strs5 = strs5 + gamma*dyy(k)*dyy(k)
        strs6 = strs6 + gamma*dyy(k)*dzz(k)
        strs9 = strs9 + gamma*dzz(k)*dzz(k)
#endif
        
c     
c     improved atomic force
        
        xxt(i)=xxt(i)+dxx(k)*gamma
        yyt(i)=yyt(i)+dyy(k)*gamma
        zzt(i)=zzt(i)+dzz(k)*gamma
        
        xxt(j)=xxt(j)-dxx(k)*gamma
        yyt(j)=yyt(j)-dyy(k)*gamma
        zzt(j)=zzt(j)-dzz(k)*gamma
        
      enddo
      
c     
c     transport temporary positions to other nodes
      
      if(mxnode.gt.1)then
        
        if(lshmov) call shmove
     x    (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x    txx,tyy,tzz,buffer)
        
      endif
      
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
        
        dli = 1.d0/dble(listme(i))
        dlj = 1.d0/dble(listme(j))

        fxx(i)=fxx(i)+xxt(i)*dli
        fyy(i)=fyy(i)+yyt(i)*dli
        fzz(i)=fzz(i)+zzt(i)*dli
        fxx(j)=fxx(j)+xxt(j)*dlj
        fyy(j)=fyy(j)+yyt(j)*dlj
        fzz(j)=fzz(j)+zzt(j)*dlj
        
      enddo
c     
c     splice force arrays across nodes
      
      if(mxnode.gt.1)then

        call gdsum(vircon,1,buffer)
        call splice 
     x    (idnode,mxnode,natms,listme,listot,fxx,fyy,fzz,buffer)
        
      endif
#ifdef STRESS
c     
c     complete stress tensor
      
      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9
#endif
  100 continue
#ifdef VAMPIR
      call VTEND(128, ierr)
#endif
      return
      end
