      subroutine pmf_shake
     x  (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x  prmpmf,virpmf,npmf,lstpmf,listpm,lstpmt,dxt,
     x  dyt,dzt,cell,pmfwght,pmfnrm,xxx,yyy,zzz,xxt,
     x  yyt,zzt,stress,summas,buffer,dxp,dyp,dzp,dsq,
     x  xa,ya,za)

c***********************************************************************
c     
c     dlpoly constraint subroutine for potential of mean force calc.
c     accummulates constraint force to maintain reaction coordinate
c     
c     copyright daresbury laboratory 1995
c     author t.forester august 1995
c     
c     wl
c     2001/05/30 12:40:21
c     1.5
c     Exp
c
c***********************************************************************

#include "dl_params.inc"

      logical safep

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension xa(2,mspmf),ya(2,mspmf),za(2,mspmf),cell(9),stress(9)
      dimension lstpmf(mxspmf,mspmf),npmf(2)
      dimension pmfwght(mxspmf),amt(2),pmfnrm(2)
      dimension listpm(mxatms),lstpmt(mxatms)
      dimension dxp(mspmf),dyp(mspmf),dzp(mspmf),dsq(mspmf)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension buffer(mxbuff),summas(2)

#ifdef VAMPIR
      call VTBEGIN(156, ierr)
#endif
      virpmf = 0.d0
      if(mxcons.lt.nspmf) call error(idnode,492)
      if(mspmf .lt.nspmf) call error(idnode,458)

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
c     assume bond vectors dxp,dyp,dzp are input
c     dxp = (sum) wght*xxx(i,1) - (sum) wght*xxx(j,2) etc

c     
c     application of constraint (shake) algorithm
      
      do icyc= 1,mxshak
c     
c     calculate bond constraint length

        do k  = 1,nspmf

          jj = 0
          kk=0 

          do ipmf = 1,2
c     
c     correct for periodic images - assume less than half box length

            i1 = lstpmf(jj+1,k)
            do i = 1,npmf(ipmf)

              jj = jj+1
              i2 = lstpmf(jj,k)
              xxt(i) = xxx(i2) - xxx(i1)
              yyt(i) = yyy(i2) - yyy(i1)
              zzt(i) = zzz(i2) - zzz(i1)
              
            enddo

            call images(imcon,0,1,npmf(ipmf),cell,xxt,yyt,zzt)
c     
c     create weighted coordinate

            xa(ipmf,k)= 0.d0
            ya(ipmf,k)= 0.d0
            za(ipmf,k)= 0.d0

            do i = 1,npmf(ipmf)

              kk = kk+1
              xa(ipmf,k) = xa(ipmf,k) + pmfwght(kk)*xxt(i)
              ya(ipmf,k) = ya(ipmf,k) + pmfwght(kk)*yyt(i)
              za(ipmf,k) = za(ipmf,k) + pmfwght(kk)*zzt(i)

            enddo

            xa(ipmf,k) = xa(ipmf,k)/pmfnrm(ipmf) + xxx(i1)
            ya(ipmf,k) = ya(ipmf,k)/pmfnrm(ipmf) + yyy(i1)
            za(ipmf,k) = za(ipmf,k)/pmfnrm(ipmf) + zzz(i1)

          enddo

          dxt(k) = xa(2,k) - xa(1,k)
          dyt(k) = ya(2,k) - ya(1,k)
          dzt(k) = za(2,k) - za(1,k)

        enddo

        call images(imcon,0,1,k,cell,dxt,dyt,dzt)

        do ipmf = 1,2
          amt(ipmf) = tstep2/summas(ipmf)*dble(1-2*mod(ipmf+1,2))
        enddo

        dis = prmpmf
        omega2 = dis*dis
        safep=.true.
        eps = 0.d0

        do k = 1,nspmf

          dsq(k) = dxt(k)**2+dyt(k)**2+dzt(k)**2
          eps = max(eps,abs((omega2-dsq(k))/(2.0d0*dis)))

        enddo
        
        if(eps.gt.tolnce) safep=.false.
        if(safep) goto 100

        do k = 1,nspmf

          gamma=(omega2-dsq(k))/(-2.d0*(amt(2)-amt(1))*
     x      (dxp(k)*dxt(k)+dyp(k)*dyt(k)+dzp(k)*dzt(k)))
c     
c     accumulate pmf virial

          virpmf=virpmf+gamma*(dxp(k)**2+dyp(k)**2+dzp(k)**2)

#ifdef STRESS
          strs1 = strs1 - gamma*dxp(k)*dxp(k)
          strs2 = strs2 - gamma*dxp(k)*dyp(k)
          strs3 = strs3 - gamma*dxp(k)*dzp(k)
          strs5 = strs5 - gamma*dyp(k)*dyp(k)
          strs6 = strs6 - gamma*dyp(k)*dzp(k)
          strs9 = strs9 - gamma*dzp(k)*dzp(k)
#endif

c     
c     improve approximate atomic positions

          jj = 0
          do ipmf=1,2

            gammi=-gamma*amt(ipmf)

            do i1 = 1,npmf(ipmf)

              jj = jj+1
              i = lstpmf(jj,k)

              xxx(i)=xxx(i)+dxp(k)*gammi
              yyy(i)=yyy(i)+dyp(k)*gammi
              zzz(i)=zzz(i)+dzp(k)*gammi

            enddo
            
          enddo

        enddo

      enddo
c
c     convergence failure - safep stays .false.

  100 continue
      
c     
c     splice coordinate arrays across nodes
      
      if(mxnode.gt.1)then
        
        call gdsum(virpmf,1,buffer)
        
        call splice 
     x    (idnode,mxnode,natms,listpm,lstpmt,xxx,yyy,zzz,buffer)
        
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
      call VTEND(156, ierr)
#endif
      return
      end
