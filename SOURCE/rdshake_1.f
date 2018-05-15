      subroutine rdshake_1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x  tolnce,tstep,vircon,lashap,lishap,listcon,listme,
     x  listot,lstfrz,buffer,cell,dxt,dxx,dyt,dyy,dzt,dzz,
     x  prmcon,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x  stress)
c     
c***********************************************************************
c     
c     dl_poly subroutine for applying bond constraint corrections after
c     atomic integration.
c     Must be used in conjunction with integration algorithms
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith august 1992.
c     t. forester march 1994.
c     
c     wl
c     2000/01/18 14:05:54
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe,lshmov
      
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
#ifdef VAMPIR
      call VTBEGIN(126, ierr)
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

c     
c     test size of work arrays

      safe =.true.
      if(mxxdf.lt.nscons) safe = .false.

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,412)

c     
c     application of constraint (shake) algorithm
      
      safe=.false.

      do icyc=1,mxshak
        
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
          
          
          dis=prmcon(listcon(k,1))

          dx =dxt(k)
          dy =dyt(k)
          dz =dzt(k)

          esig1=(abs(dx*dx + dy*dy + dz*dz - dis*dis)/
     x      (dis))
          esig = max(esig,esig1)

        enddo
        
        esig = esig*0.5d0
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
c     initialise increment arrays

        do i=1,natms

          xxt(i)=0.d0
          yyt(i)=0.d0
          zzt(i)=0.d0

        enddo

c     
c     calculate constraint forces
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
c     
c     set constraint parameters
          
          dis=prmcon(listcon(k,1))
          omega2= dis*dis
          amti= tstep2/weight(i)
          amtj=-tstep2/weight(j)
          
          if(lstfrz(i).ne.0) amti = 0.d0
          if(lstfrz(j).ne.0) amtj = 0.d0
c     
c     constraint force parameter
          
          dx = dxt(k)
          dy = dyt(k)
          dz = dzt(k)
          
          gamma=(omega2-(dx*dx+dy*dy+dz*dz))/
     x      (-2.d0*(amti-amtj)*
     x      (dxx(k)*dx+dyy(k)*dy+dzz(k)*dz))

c     
c     accumulate bond virial

          vircon=vircon+gamma*(dxx(k)**2+dyy(k)**2+dzz(k)**2)
#ifdef STRESS
          strs1 = strs1 - gamma*dxx(k)*dxx(k)
          strs2 = strs2 - gamma*dxx(k)*dyy(k)
          strs3 = strs3 - gamma*dxx(k)*dzz(k)
          strs5 = strs5 - gamma*dyy(k)*dyy(k)
          strs6 = strs6 - gamma*dyy(k)*dzz(k)
          strs9 = strs9 - gamma*dzz(k)*dzz(k)
#endif
c     
c     improve approximate atomic positions
          
          gammi=-gamma*amti
          xxt(i)=xxt(i)+dxx(k)*gammi
          yyt(i)=yyt(i)+dyy(k)*gammi
          zzt(i)=zzt(i)+dzz(k)*gammi

          gammj=-gamma*amtj
          xxt(j)=xxt(j)+dxx(k)*gammj
          yyt(j)=yyt(j)+dyy(k)*gammj
          zzt(j)=zzt(j)+dzz(k)*gammj
          
        enddo
        
c     
c     transport temporary positions to other nodes
        
        if(mxnode.gt.1)then
          
          if(lshmov) call shmove
     x      (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x      txx,tyy,tzz,buffer)

        endif
        
        do k=1,nscons

          i=listcon(k,2)
          j=listcon(k,3)
          
          dli = 1.d0/dble(listme(i))
          dlj = 1.d0/dble(listme(j))

          xxx(i)=xxx(i)+xxt(i)*dli
          yyy(i)=yyy(i)+yyt(i)*dli
          zzz(i)=zzz(i)+zzt(i)*dli
          xxx(j)=xxx(j)+xxt(j)*dlj
          yyy(j)=yyy(j)+yyt(j)*dlj
          zzz(j)=zzz(j)+zzt(j)*dlj

        enddo

      enddo
      
c     
c     error exit for non-convergence

      return
      
  100 continue

c     
c     splice coordinate arrays across nodes

      if(mxnode.gt.1)then

        call gdsum(vircon,1,buffer)
        call splice 
     x    (idnode,mxnode,natms,listme,listot,xxx,yyy,zzz,buffer)

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
#ifdef VAMPIR
      call VTEND(126, ierr)
#endif
      return
      end
