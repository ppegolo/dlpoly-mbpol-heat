      subroutine pmf_1
     x  (safe,safep,lshmov,idnode,imcon,mxnode,natms,nscons,
     x  ntcons,nspmf,ntpmf,engke,prmpmf,tolnce,tstep,vircon,
     x  virpmf,lstpmf,npmf,lashap,lishap,listcon,listpm,lstpmt,
     x  listme,listot,lstfrz,buffer,cell,dxt,dxx,dyt,
     x  dyy,dzt,dzz,pmfwght,fxx,fyy,fzz,prmcon,txx,tyy,tzz,
     x  uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,xxx,ydf,yyt,
     x  yyy,zdf,zzt,zzz,stress,xa,ya,za,dxp,dyp,dzp,dsq)
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics. Verlet leapfrog With RD-SHAKE
c     and PMF_SHAKE - for potential of mean force calculations.
c     
c     parallel replicated data version : block data
c     adapted from dl_poly routine nve_1.f
c     
c     copyright - daresbury laboratory 1995
c     author  - t.forester aug 1995
c     
c     wl
c     2001/05/30 12:40:21
c     1.7
c     $Sate: Exp $
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe,lshmov,safep
      
      dimension listme(mxatms),lishap(mxlshp),lstfrz(mxatms)
      dimension lashap(mxproc),listot(mxatms)
      dimension listcon(mxcons,3),buffer(mxbuff)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension weight(mxatms),prmcon(mxtcon),cell(9)
      dimension stress(9),stres1(9),dsq(mspmf)
      dimension lstpmf(mxspmf,mspmf),npmf(2),pmfnrm(2)
      dimension pmfwght(mxspmf),summas(2)
      dimension xa(2,mspmf),ya(2,mspmf),za(2,mspmf)
      dimension dxp(mspmf),dyp(mspmf),dzp(mspmf)
      dimension listpm(mxpmf),lstpmt(mxpmf)
      
#ifdef VAMPIR
      call VTBEGIN(58, ierr)
#endif
      nbuff=mxbuff
c     
c     block indices

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
c     
c     constraint virials

      vircon=0.d0
      virpmf=0.d0

#ifdef STRESS
c     
c     temporary stress tensor accumulators
      do i = 1,9
        stres1(i) = 0.d0
      enddo
#endif
c     
c     store initial values of position
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xdf(j)=xxx(i)
        ydf(j)=yyy(i)
        zdf(j)=zzz(i)
        
      enddo
c     
c     construct current bond vectors
      
      do k=1,nscons
c     
c     indices of atoms in bond
        
        i=listcon(k,2)
        j=listcon(k,3)
c     
c     calculate current bond vector
        
        dxx(k)=xxx(i)-xxx(j)
        dyy(k)=yyy(i)-yyy(j)
        dzz(k)=zzz(i)-zzz(j)
        
      enddo
c     
c     periodic boundary condition for bond vectors
      
      call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
c     
c     calculate PMF bond constraints and store initial positions

      do k = 1,nspmf

        jj = 0
        kk = 0
        
        do ipmf = 1,2
c     
c     correct for periodic images - assume less than half box length

          i1 = lstpmf(jj+1,k)
          summas(ipmf) = 0.d0
          do i = 1,npmf(ipmf)

            jj = jj+1
            i2 = lstpmf(jj,k)
            xxt(i) = xxx(i2) - xxx(i1)
            yyt(i) = yyy(i2) - yyy(i1)
            zzt(i) = zzz(i2) - zzz(i1)
            summas(ipmf) = summas(ipmf) + weight(i2)

          enddo

          call images(imcon,0,1,npmf(ipmf),cell,xxt,yyt,zzt)
c     
c     create weighted coordinate

          xa(ipmf,k)= 0.d0
          ya(ipmf,k)= 0.d0
          za(ipmf,k)= 0.d0
          pmfnrm(ipmf)= 0.d0
          
          do i = 1,npmf(ipmf)

            kk = kk+1
            xa(ipmf,k) = xa(ipmf,k) + pmfwght(kk)*xxt(i)
            ya(ipmf,k) = ya(ipmf,k) + pmfwght(kk)*yyt(i)
            za(ipmf,k) = za(ipmf,k) + pmfwght(kk)*zzt(i)
            pmfnrm(ipmf) = pmfnrm(ipmf) + pmfwght(kk)

          enddo

          xa(ipmf,k) = xa(ipmf,k)/pmfnrm(ipmf) + xxx(i1)
          ya(ipmf,k) = ya(ipmf,k)/pmfnrm(ipmf) + yyy(i1)
          za(ipmf,k) = za(ipmf,k)/pmfnrm(ipmf) + zzz(i1)

        enddo


        dxp(k) = xa(2,k) - xa(1,k)
        dyp(k) = ya(2,k) - ya(1,k)
        dzp(k) = za(2,k) - za(1,k)

      enddo
c     
c     periodic boundary condition for pmf bond vectors
      
      call images(imcon,0,1,nspmf,cell,dxp,dyp,dzp)

c
c     move atoms by leapfrog algorithm

      safe=(ntcons.eq.0)
      safep=(ntpmf.eq.0)
      j = 0
      do i=iatm0,iatm1
        j = j+1
c     
c     update velocities
        
        uxx(i)=vxx(i) + tstep*fxx(i)/weight(i)
        uyy(i)=vyy(i) + tstep*fyy(i)/weight(i)
        uzz(i)=vzz(i) + tstep*fzz(i)/weight(i)
c     
c     update positions
        
        xxx(i)=xdf(j)+tstep*uxx(i)
        yyy(i)=ydf(j)+tstep*uyy(i)
        zzz(i)=zdf(j)+tstep*uzz(i)
        
      enddo
      
      if(ntcons.gt.0.or.ntpmf.gt.0) then
c     
c     RDSHAKE procedure 
c     
c     global exchange of configuration data
        if(mxnode.gt.1) then
          
          call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
          
        endif

c     
c     apply constraint corrections - iteratively

        do icyc = 1,mxshak
c
c     apply bond constraints

          viracc = 0.d0
          if(ntcons.gt.0) call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,lashap,lishap,listcon,listme,
     x      listot,lstfrz,buffer,cell,dxt,dxx,dyt,dyy,dzt,dzz,
     x      prmcon,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x      stres1)

          vircon = vircon+viracc
c     
c     apply pmf constraints

          virac1 = 0.d0

          if(ntpmf.gt.0) call pmf_shake
     x      (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x      prmpmf,virac1,npmf,lstpmf,listpm,lstpmt,dxt,
     x      dyt,dzt,cell,pmfwght,pmfnrm,xxx,yyy,zzz,xxt,
     x      yyt,zzt,stres1,summas,buffer,dxp,dyp,dzp,dsq,
     x      xa,ya,za)

          virpmf = virpmf + virac1
c
c     exit loop if converged

          if(safe.and.safep.and.virac1.eq.0.d0)goto 100
          
        enddo

        safep=.false.

  100   continue
c     
c     calculate velocity correction
        
        rstep = 1.d0/tstep
        j=0
        do i=iatm0,iatm1
          
c     
c     update corrected velocity
          
          j=j+1
          uxx(i)=(xxx(i)-xdf(j))*rstep
          uyy(i)=(yyy(i)-ydf(j))*rstep
          uzz(i)=(zzz(i)-zdf(j))*rstep
          
c     
c     calculate the corrected forces
          
          fxx(i)=(uxx(i)-vxx(i))*weight(i)*rstep
          fyy(i)=(uyy(i)-vyy(i))*weight(i)*rstep
          fzz(i)=(uzz(i)-vzz(i))*weight(i)*rstep
          
        enddo
        

      endif
c     
c     updated velocity and kinetic energy
      
      engke =0.d0
      do i = iatm0,iatm1

        vxt = 0.5d0*(vxx(i)+uxx(i))
        vyt = 0.5d0*(vyy(i)+uyy(i))
        vzt = 0.5d0*(vzz(i)+uzz(i))

        engke=engke+0.5d0*weight(i)*(vxt*vxt+vyt*vyt+vzt*vzt)
        
        vxx(i) = uxx(i)
        vyy(i) = uyy(i)
        vzz(i) = uzz(i)
#ifdef STRESS
c     
c     kinetic contribution to stress tensor

        stres1(1) = stres1(1) + weight(i)*vxt*vxt
        stres1(2) = stres1(2) + weight(i)*vxt*vyt
        stres1(3) = stres1(3) + weight(i)*vxt*vzt
        stres1(5) = stres1(5) + weight(i)*vyt*vyt
        stres1(6) = stres1(6) + weight(i)*vyt*vzt
        stres1(9) = stres1(9) + weight(i)*vzt*vzt
#endif
      enddo
c     
c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c     
c     global exchange of configuration data
      
      if(mxnode.gt.1)then

        nbuff=mxbuff
        call gdsum(engke,1,buffer)
        call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,nbuff,vxx,vyy,vzz,buffer)
        
        if(ntcons.gt.0) then
          
          call merge(idnode,mxnode,natms,nbuff,fxx,fyy,fzz,buffer)
          
        endif
        
      endif
#ifdef STRESS
c     
c     complete stress tensor

      if(mxnode.gt.1) call gdsum(stres1,9,buffer)

      stress(1) = stress(1) + stres1(1)
      stress(2) = stress(2) + stres1(2)
      stress(3) = stress(3) + stres1(3)
      stress(4) = stress(4) + stres1(2)
      stress(5) = stress(5) + stres1(5)
      stress(6) = stress(6) + stres1(6)
      stress(7) = stress(7) + stres1(3)
      stress(8) = stress(8) + stres1(6)
      stress(9) = stress(9) + stres1(9)
#endif
#ifdef VAMPIR
      call VTEND(58, ierr)
#endif
      return
      end



