      subroutine nvt_b0
     x  (idnode,mxnode,imcon,natms,engke,sigma,tstep,taut,
     x  buffer,cell,fxx,fyy,fzz,vxx,vyy,vzz,weight,xxx,yyy,
     x  zzz,stress)
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory  1993
c     author    - t. forester      july 1993
c     
c     wl
c     2000/01/18 14:05:48
c     1.3
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
      
c     
c     set kinetic energy accumulator
      
#ifdef VAMPIR
      call VTBEGIN(48, ierr)
#endif
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
c     estimate current temperature at half step
        
        vxt = vxx(i) + tstep*fxx(i)/weight(i)*0.5d0
        vyt = vyy(i) + tstep*fyy(i)/weight(i)*0.5d0
        vzt = vzz(i) + tstep*fzz(i)/weight(i)*0.5d0
        
        engke = engke + weight(i)*(vxt*vxt + vyt*vyt + vzt*vzt)
        
      enddo
      
      if(mxnode.gt.1) call gdsum(engke,1,buffer)
      
      engke = engke*0.5d0
      
c     
c     temperature scaling  coefficient - taut is the decay constant
      
      chi =  sqrt(1.d0 + tstep/taut*(sigma/engke -1.d0))
      
      engke = 0.0d0
      
c     move atoms by leapfrog algorithm
      
      do i=iatm0,iatm1
        
        uxx = vxx(i)
        uyy = vyy(i)
        uzz = vzz(i)
c     
c     calculate new velocities
        
        vxx(i)=chi*(vxx(i) + tstep*(fxx(i)/weight(i)))
        vyy(i)=chi*(vyy(i) + tstep*(fyy(i)/weight(i)))
        vzz(i)=chi*(vzz(i) + tstep*(fzz(i)/weight(i)))
c
c     velocity at half time step

        vxt = 0.5d0*(vxx(i)+uxx)
        vyt = 0.5d0*(vyy(i)+uyy)
        vzt = 0.5d0*(vzz(i)+uzz)
c     
c     update positions
        
        xxx(i)=xxx(i)+tstep*vxx(i)
        yyy(i)=yyy(i)+tstep*vyy(i)
        zzz(i)=zzz(i)+tstep*vzz(i)
c     
c     calculate kinetic energy*2
        
        engke=engke+weight(i)*(vxt*vxt+vyt*vyt+vzt*vzt)

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
c     kinetic energy
      
      engke = engke*0.5d0
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
      call VTEND(48, ierr)
#endif
      return
      end
