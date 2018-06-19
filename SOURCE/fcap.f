      subroutine fcap
     x  (lfcap,idnode,mxnode,natms,fmax,temp,fxx,fyy,fzz)
      
c*********************************************************************
c     
c     DLPOLY routinue for limiting the absolute magnitude of
c     forces. Used in equilibration period only
c     
c     copyright daresbury laboratory 1993
c     
c     author -     t. forester march 1993
c     amended-     t. forester  sept 1994
c     
c     wl
c     2000/01/18 14:05:39
c     1.3
c     Exp
c     
c*********************************************************************
c     
#include "dl_params.inc"
      
      logical lfcap
      
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      
#ifdef VAMPIR
      call VTBEGIN(29, ierr)
#endif
      if(lfcap) then
c     
c     maximum force permitted
        
        fmax1 = boltz*fmax*temp
        fmax2 = fmax1*fmax1
c     
c     cap forces and conserve linear momentum
        
        fxc = 0.d0
        fyc = 0.d0
        fzc = 0.d0
        
        do i = 1,natms
          
          fmod = fxx(i)**2 + fyy(i)**2 + fzz(i)**2
          
          if(fmod.gt.fmax2) then
            
            fscale = sqrt(fmax2/fmod)
            
            fxx(i) = fxx(i)*fscale
            fyy(i) = fyy(i)*fscale
            fzz(i) = fzz(i)*fscale
            
          endif
c     
c     accummulate forces - to check on momentum conservation
          
          fxc = fxc + fxx(i)
          fyc = fyc + fyy(i)
          fzc = fzc + fzz(i)
          
        enddo
c     
c     ensure net forces sum to zero
        
        fxc = -fxc/dble(natms)
        fyc = -fyc/dble(natms)
        fzc = -fzc/dble(natms)
c     
c     conserve momentum
        
        do i = 1,natms
          
          fxx(i) = fxx(i) + fxc
          fyy(i) = fyy(i) + fyc
          fzz(i) = fzz(i) + fzc
          
        enddo
        
      endif
#ifdef VAMPIR
      call VTEND(29, ierr)
#endif
      return
      end
