      subroutine nve_0
     x  (idnode,mxnode,natms,imcon,tstep,engke,weight,cell,
     x  xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,buffer,stress)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     
c     wl
c     2000/01/18 14:05:47
c     1.4
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
      call VTBEGIN(32, ierr)
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
c     block indices
      
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode
      
c     
c     move atoms by leapfrog algorithm
      
      do i=iatm0,iatm1
c     
c     store old velocities
        
        uxx=vxx(i)
        uyy=vyy(i)
        uzz=vzz(i)
c     
c     calculate new velocities
        
        vxx(i)=vxx(i)+(tstep/weight(i))*fxx(i)
        vyy(i)=vyy(i)+(tstep/weight(i))*fyy(i)
        vzz(i)=vzz(i)+(tstep/weight(i))*fzz(i)
        
c     
c     update positions
        
        xxx(i)=xxx(i)+tstep*vxx(i)
        yyy(i)=yyy(i)+tstep*vyy(i)
        zzz(i)=zzz(i)+tstep*vzz(i)
        
c     
c     calculate kinetic energy
        
        engke=engke+(weight(i)/8.d0)*
     x    ((vxx(i)+uxx)**2+(vyy(i)+uyy)**2+(vzz(i)+uzz)**2)
        
#ifdef STRESS
c     
c     kinetic contribution to stress tensor
        
        vxt=0.5d0*(vxx(i)+uxx)
        vyt=0.5d0*(vyy(i)+uyy)
        vzt=0.5d0*(vzz(i)+uzz)
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
      call VTEND(32, ierr)
#endif
      return
      end
