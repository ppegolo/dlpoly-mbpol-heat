      subroutine freeze
     x  (idnode,mxnode,natms,lstfrz,vxx,vyy,vzz,fxx,fyy,fzz)

c***********************************************************************
c     
c     dlpoly routine to quench forces and velocities on 'frozen' atoms
c     replicated data version - blocked data
c     
c     copyright daresbury laboratory 1994
c     author t.forester nov 1994
c     
c     wl
c     2000/01/18 14:05:39
c     1.4
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      dimension lstfrz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
#ifdef VAMPIR
      call VTBEGIN(31, ierr)
#endif
c     
c     block indices

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode
      
      do i = iatm1,iatm2
        
        if(lstfrz(i).ne.0) then
          
          vxx(i) = 0.d0
          vyy(i) = 0.d0
          vzz(i) = 0.d0
          fxx(i) = 0.d0
          fyy(i) = 0.d0
          fzz(i) = 0.d0
          
        endif
        
      enddo

#ifdef VAMPIR
      call VTEND(31, ierr)
#endif
      return
      end
