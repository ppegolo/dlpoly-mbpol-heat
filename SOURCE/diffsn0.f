      subroutine diffsn0
     x  (idnode,natms,mxnode,tstep,vxx,vyy,vzz,xx0,yy0,zz0)

c***********************************************************************
c     
c     DL_POLY routine for calculating displacements of sites from
c     t=0 positions
c     
c     use diffsn1 for mean squared displacements
c     
c     parallel version - replicated data.
c     
c     copyright daresbury laboratory 1993
c     
c     author - t. forester      june 1993
c     
c     wl
c     2000/01/18 14:05:34
c     1.3
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      logical newjob 

      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension xx0(mxatms),yy0(mxatms),zz0(mxatms)

      save newjob,iatm1,iatm2
      data newjob/.true./

#ifdef VAMPIR
      call VTBEGIN(74, ierr)
#endif
      if(newjob) then

        newjob=.false.
        iatm1 = (idnode*natms)/mxnode+1
        iatm2 = ((idnode+1)*natms)/mxnode

      endif

      do i = iatm1,iatm2

        xx0(i) = xx0(i) + vxx(i)*tstep
        yy0(i) = yy0(i) + vyy(i)*tstep
        zz0(i) = zz0(i) + vzz(i)*tstep

      enddo

#ifdef VAMPIR
      call VTEND(74, ierr)
#endif
      return
      end

