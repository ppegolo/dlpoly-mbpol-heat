      subroutine nvt_e1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x  engke,sigma,tolnce,tstep,vircon,lashap,lishap,
     x  listcon,listme,listot,lstfrz,buffer,cell,dxt,
     x  dxx,dyt,dyy,dzt,dzz,fxx,fyy,fzz,prmcon,txx,tyy,
     x  tzz,uxx,uyy,uzz,vxx,vyy,vzz,weight,xdf,xxt,xxx,
     x  ydf,yyt,yyy,zdf,zzt,zzz,stress)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Evans
c     thermostat.
c     Comp. Phys. reports 1, 299, (1984)
c     
c     parallel replicated data version : block data
c     
c     for systems using bond CONSTRAINTS.
c     
c     copyright - daresbury laboratory 1993
c     author    -    t. forester   july 1993
c     amended   -    t. forester   sept 1994
c     amended   -    t. forester   dec  1994 : block data
c     
c     wl
c     2000/01/18 14:05:49
c     1.4
c     Exp
c
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical safe,lshmov
      
      dimension listme(mxatms),lishap(mxlshp),lashap(mxproc)
      dimension listot(mxatms),listcon(mxcons,3),lstfrz(mxatms)
      dimension buffer(mxbuff)
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
      dimension stress(9),stres1(9)
      
#ifdef VAMPIR
      call VTBEGIN(51, ierr)
#endif
      safe=.false.
      nbuff=mxbuff
c     
c     block indices

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode
c     
c     constarint virial

      vircon=0.d0
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
      do i=iatm1,iatm2
        
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
c     estimate kinetic energy: at half step
      
      engke = 0.d0
      j = 0
      
      do i = iatm1,iatm2
        
        j = j+1
c     
c     estimate new velocity at the half step
        
        tmp = tstep*0.5d0/weight(i)
        vxt = vxx(i) + tmp*fxx(i)
        vyt = vyy(i) + tmp*fyy(i)
        vzt = vzz(i) + tmp*fzz(i)
        
        engke = engke + weight(i)*(vxt*vxt+vyt*vyt+vzt*vzt)
        
      enddo
      
      engke = engke*0.5d0
      
      if(mxnode.gt.1) call gdsum(engke,1,buffer)
      
c     
c     begin temperature control iteration
      
      mxiter = 2
      if(ntcons.eq.0) mxiter = 1
      
      do iter = 1,mxiter
c     
c     friction coefficient for thermostat
        
        chi = sqrt(sigma/engke)
        twochi = 2.0d0*chi - 1.0d0
c     
c     move atoms by leapfrog algorithm
        
        engke = 0.d0
        j = 0
        do i=iatm1,iatm2
          j = j+1
c     
c     update velocities
          
          tchi = chi*tstep/weight(i)         
          uxx(i)=vxx(i)*twochi + tchi*fxx(i)
          uyy(i)=vyy(i)*twochi + tchi*fyy(i)
          uzz(i)=vzz(i)*twochi + tchi*fzz(i)
c     
c     update positions
          
          xxx(i)=xdf(j)+tstep*uxx(i)
          yyy(i)=ydf(j)+tstep*uyy(i)
          zzz(i)=zdf(j)+tstep*uzz(i)
          
        enddo
        
        if(ntcons.eq.0) safe =.true.
        if(ntcons.gt.0) then
c     
c     Begin shake procedures
c     
c     global exchange of configuration data
          
          if(mxnode.gt.1) then
            
            call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
            
          endif
          
c     
c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,lashap,lishap,listcon,listme,
     x      listot,lstfrz,buffer,cell,dxt,dxx,dyt,dyy,dzt,dzz,
     x      prmcon,txx,tyy,tzz,weight,xxt,xxx,yyt,yyy,zzt,zzz,
     x      stres1)
c     
c     contribution to constraint virial 
          
          vircon = vircon+viracc
c     
c     calculate velocity correction
          
          rstep = 1.d0/tstep
          j=0
          do i=iatm1,iatm2
            
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
c     kinetic energy
        
        engke =0.d0
        
        do i=iatm1,iatm2
          
          engke = engke+weight(i)*((vxx(i)+uxx(i))**2
     x      + (vyy(i)+uyy(i))**2 +(vzz(i)+uzz(i))**2)
          
          
        enddo
        
        engke = engke*0.125d0
        if(mxnode.gt.1) call gdsum(engke,1,buffer)
c     
c     end of thermal constraint iteration
        
      enddo
      
      engke = 0.d0
c     
c     updated velocity

      do i = iatm1,iatm2

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
      call VTEND(51, ierr)
#endif
      return
      end
