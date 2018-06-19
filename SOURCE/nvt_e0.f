      subroutine nvt_e0
     x  (idnode,imcon,natms,mxnode,engke,sigma,tstep,buffer,
     x  cell,fxx,fyy,fzz,vxx,vyy,vzz,weight,xxx,yyy,zzz,stress)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Evans
c     thermostat.
c     Comp. Phys. reports 1, 299, (1984)
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - t. forester   feb 1993
c
c     wl
c     2000/01/18 14:05:49
c     1.5
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension weight(mxatms),buffer(mxbuff),cell(9)
      dimension stress(9),stres1(9)

#ifdef VAMPIR
      call VTBEGIN(50, ierr)
#endif
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
        
        vxt = vxx(i) + tstep*fxx(i)/weight(i)*0.5d0
        vyt = vyy(i) + tstep*fyy(i)/weight(i)*0.5d0
        vzt = vzz(i) + tstep*fzz(i)/weight(i)*0.5d0
        
        engke = engke + weight(i)*(vxt**2+
     $    vyt**2 + vzt**2)
        
      enddo
      
      engke = engke*0.5d0
      
      if(mxnode.gt.1) call gdsum(engke,1,buffer)
      
      chi = sqrt(sigma/engke)
      twochi = 2.0d0*chi - 1.0d0
      chi = chi*tstep
      
c     move atoms by leapfrog algorithm
      
      engke = 0.0d0
      
      do i=iatm0,iatm1
c     
c     calculate new velocities
        
        uxx = vxx(i)
        uyy = vyy(i)
        uzz = vzz(i)
c     
c     update velocities
        
        tchi = chi/weight(i)         
        vxx(i)=vxx(i)*twochi + tchi*fxx(i)
        vyy(i)=vyy(i)*twochi + tchi*fyy(i)
        vzz(i)=vzz(i)*twochi + tchi*fzz(i)
c     
c     update positions
        
        xxx(i)=xxx(i)+tstep*vxx(i)
        yyy(i)=yyy(i)+tstep*vyy(i)
        zzz(i)=zzz(i)+tstep*vzz(i)
        
c     
c     calculate kinetic energy

        vxt = 0.5d0*(vxx(i)+uxx)
        vyt = 0.5d0*(vyy(i)+uyy)
        vzt = 0.5d0*(vzz(i)+uzz)

        engke=engke+0.5d0*weight(i)*(vxt*vxt+vyt*vyt+vzt*vzt)
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
      call VTEND(50, ierr)
#endif
      return
      end

