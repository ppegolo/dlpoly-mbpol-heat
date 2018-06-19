      subroutine nvt_h0
     x  (idnode,mxnode,imcon,natms,chit,conint,consv,engke,
     x  taut,sigma,tstep,buffer,cell,fxx,fyy,fzz,vxx,vyy,
     x  vzz,weight,xxx,yyy,zzz,uxx,uyy,uzz,stress)
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover
c     thermostat.
c     Phys. Rev A 31, 1695 (1985)
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester      feb 1993
c     amended   - t. forester     sept 1994
c     
c     wl
c     2001/05/30 12:40:19
c     1.6
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension weight(mxatms),buffer(mxbuff),cell(9)
      dimension stress(9),stres1(9)
#ifdef VAMPIR
      call VTBEGIN(52, ierr)
#endif
c     
c     inertia parameter for Nose Hoover thermostat
      
      qmass = 2.0d0*sigma*taut**2
      
c     
c     set kinetic energy accumulator
      
      engke=0.d0
#ifdef STRESS
c     
c     temporary stress tensor accumulators
      do i = 1,9
        stres1(i) = 0.d0
      enddo
#endif
c     
c     set up block indices
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
      
      do i = iatm0,iatm1
c     
c     estimate current temperature
        
        uxx(i) = vxx(i) + tstep*fxx(i)/weight(i)*0.5d0
        uyy(i) = vyy(i) + tstep*fyy(i)/weight(i)*0.5d0
        uzz(i) = vzz(i) + tstep*fzz(i)/weight(i)*0.5d0
        
        engke = engke + weight(i)*(uxx(i)**2+
     x    uyy(i)**2 + uzz(i)**2)
        
      enddo
      
      engke = engke*0.5d0
      
      if(mxnode.gt.1) call gdsum(engke,1,buffer)
c     
c     begin thermostat iteration
      
      maxit = 3
      do iter = 1,maxit
        
c     
c     time derivative of friction coefficient
        
        chip = 2.d0*(engke - sigma)/qmass
        
c     
c     propagate friction coefficient
        
        chinew =  chit + tstep*chip
c     
c     value at half time step
        
        chit1 = 0.5d0*(chinew+chit)
c     
c     kinetic energy accumulator
        
        engke = 0.0d0
        
        if(iter.ne.maxit) then
c     
c     find better estimate to kinetic energy at half timestep
          
          do i=iatm0,iatm1
c     
c     calculate new velocities at half step
            
            uxx(i)=vxx(i)+0.5d0*tstep*(fxx(i)/weight(i)-chit1*uxx(i))
            uyy(i)=vyy(i)+0.5d0*tstep*(fyy(i)/weight(i)-chit1*uyy(i))
            uzz(i)=vzz(i)+0.5d0*tstep*(fzz(i)/weight(i)-chit1*uzz(i))
c     
c     calculate kinetic energy
            
            engke=engke+weight(i)*(uxx(i)**2+uyy(i)**2+uzz(i)**2)
            
          enddo
          
          engke = engke*0.5d0
          
        else
c     
c     final update of velocities and positions
          
          do i = iatm0,iatm1
            
            uxx(i)=vxx(i)+tstep*(fxx(i)/weight(i)-chit1*uxx(i))
            uyy(i)=vyy(i)+tstep*(fyy(i)/weight(i)-chit1*uyy(i))
            uzz(i)=vzz(i)+tstep*(fzz(i)/weight(i)-chit1*uzz(i))
            
            vxt = 0.5d0*(uxx(i)+vxx(i))
            vyt = 0.5d0*(uyy(i)+vyy(i))
            vzt = 0.5d0*(uzz(i)+vzz(i))
c     
c     calculate kinetic energy
            
            engke=engke+weight(i)*(vxt*vxt + vyt*vyt + vzt*vzt)
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
c     
c     update velocities
            
            vxx(i) = uxx(i)
            vyy(i) = uyy(i)
            vzz(i) = uzz(i)
c     
c     update positions
            
            xxx(i)=xxx(i)+tstep*vxx(i)
            yyy(i)=yyy(i)+tstep*vyy(i)
            zzz(i)=zzz(i)+tstep*vzz(i)
            
          enddo
          
          engke = engke*0.5d0
          
        endif
c     
c     global sum of kinetic energy
        
        if(mxnode.gt.1) call gdsum(engke,1,buffer)
        
      enddo
      
c     
c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c     
c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        nbuff=mxbuff
        call merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,nbuff,vxx,vyy,vzz,buffer)
        
      endif
c     
c     update thermostat variable

      chit = chitnew
c     
c     integral of friction coefficient
      
      conint = conint + tstep*2.d0*sigma*chit
c     
c     thermostating rate squared
      
      vt2 = 2.d0*sigma/qmass
c     
c     conserved quantity less kinetic and potential energy terms
      
      consv = conint+ sigma*(chit**2/vt2)
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
      call VTEND(52, ierr)
#endif
      return
      end


