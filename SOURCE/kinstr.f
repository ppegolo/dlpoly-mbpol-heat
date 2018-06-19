      subroutine kinstr
     x  (idnode,mxnode,natms,tstep,vxx,vyy,vzz,
     x  fxx,fyy,fzz,weight,stress)
      
c***********************************************************************
c     
c     dlpoly routine to calculate the kinetic energy contribution to
c     the stress tensor
c     
c     assumes velocities are half-timestep behind forces
c     
c     replicated data version / block data
c     
c     copyright daresbury laboratory 1994
c     author t.forester may 1994
c     amended t.forester dec 1994 : block data
c     
c     wl
c     2001/05/30 12:40:09
c     1.2
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension stress(9),weight(mxatms)
#ifdef VAMPIR
      call VTBEGIN(146, ierr)
#endif
c     
c     block indices

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      do i = iatm1,iatm2

        if(weight(i).gt.0.d0) then

          vxt = vxx(i)+fxx(i)/weight(i)*tstep*0.5d0
          vyt = vyy(i)+fyy(i)/weight(i)*tstep*0.5d0
          vzt = vzz(i)+fzz(i)/weight(i)*tstep*0.5d0

          stress(1)=stress(1)-weight(i)*vxt*vxt
          stress(2)=stress(2)-weight(i)*vxt*vyt
          stress(3)=stress(3)-weight(i)*vxt*vzt

          stress(4)=stress(4)-weight(i)*vyt*vxt
          stress(5)=stress(5)-weight(i)*vyt*vyt
          stress(6)=stress(6)-weight(i)*vyt*vzt

          stress(7)=stress(7)-weight(i)*vzt*vxt
          stress(8)=stress(8)-weight(i)*vzt*vyt
          stress(9)=stress(9)-weight(i)*vzt*vzt

        endif

      enddo

#ifdef VAMPIR
      call VTEND(146, ierr)
#endif
      return
      end
