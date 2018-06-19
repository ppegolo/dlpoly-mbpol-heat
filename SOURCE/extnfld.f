      subroutine extnfld
     x  (idnode,imcon,keyfld,mxnode,natms,engfld,virfld,cell,
     x  chge,prmfld,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,weight)

c***********************************************************************
c     
c     dl_poly routine for application of an external field
c     
c     replicated data version / block data
c     
c     copyright daresbury laboratory 1993
c     author -    t.forester october 1993
c     amended-    t.forester dec 1994
c     
c     wl
c     2001/05/30 12:40:05
c     1.4
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension chge(mxatms),cell(9),weight(mxatms)
      dimension prmfld(mxfld)
#ifdef VAMPIR
      call VTBEGIN(27, ierr)
#endif
c     
c     energy and virial accumulators 
      engfld = 0.d0
      virfld = 0.d0
c     
c     block indices

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      if(keyfld.eq.1) then
c     
c     electric field: prmfld(1-3) are field components

        do i = iatm1,iatm2

          fxx(i) = fxx(i)+ chge(i)*prmfld(1)
          fyy(i) = fyy(i)+ chge(i)*prmfld(2)
          fzz(i) = fzz(i)+ chge(i)*prmfld(3)

        enddo

      elseif(keyfld.eq.2) then
c     
c     oscillating shear: orthorhombic box:  Fx = a*cos(b.2.pi.z/L)

        rz = 2.d0*pi/cell(9)

        do i = iatm1,iatm2

          fxx(i) = fxx(i) + prmfld(1)*cos(prmfld(2)*zzz(i)*rz)
          
        enddo

      elseif(keyfld.eq.3.and.imcon.eq.6) then
c     
c     continuous shear of walls : 2D periodic box (imcon=6)
c     shear rate = prmfld(1) angstrom per ps for atoms at
c     abs(z) > prmfld(2)

        do i=iatm1,iatm2

          if(abs(zzz(i)).gt.prmfld(2)) then

            vxx(i) = 0.5d0*sign(prmfld(1),zzz(i))

          endif

        enddo

      elseif(keyfld.eq.4) then
c     
c     gravitational field: field components given by prmfld(1-3)

        do i =iatm1,iatm2

          fxx(i) = fxx(i) + prmfld(1)*weight(i)
          fyy(i) = fyy(i) + prmfld(2)*weight(i)
          fzz(i) = fzz(i) + prmfld(3)*weight(i)

        enddo

      elseif(keyfld.eq.5) then

c     
c     magnetic field: field components given by prmfld(1-3)

        do i = iatm1,iatm2

          fxx(i)=fxx(i)+(vyy(i)*prmfld(3)-vzz(i)*prmfld(2))
     x      *chge(i)
          fyy(i)=fyy(i)+(vzz(i)*prmfld(1)-vxx(i)*prmfld(3))
     x      *chge(i)
          fzz(i)=fzz(i)+(vxx(i)*prmfld(2)-vyy(i)*prmfld(1))
     x      *chge(i)

        enddo

      elseif(keyfld.eq.6) then

c     
c     containing sphere : r^(-n) potential

        do i = iatm1,iatm2

          rrr = sqrt(xxx(i)**2+yyy(i)**2+zzz(i)**2)
          if(rrr.gt.prmfld(4)) then
            rrr = prmfld(2) - rrr
            if(rrr.lt.0.d0) rrr = 0.1d0

            gamma  = prmfld(1)*rrr**(-prmfld(3))
            engfld = engfld + gamma

            gamma = -prmfld(3)*gamma/((prmfld(2)-rrr)*rrr)

            fxx(i)=fxx(i)+ gamma*xxx(i)
            fyy(i)=fyy(i)+ gamma*yyy(i)
            fzz(i)=fzz(i)+ gamma*zzz(i)

          endif

        enddo

      elseif(keyfld.eq.7) then

c     
c     repulsive wall (harmonic) starting at z0

        do i = iatm1,iatm2

          if(prmfld(3)*zzz(i).gt.prmfld(2)) then

            zdif = zzz(i) - prmfld(2)
            gamma = -prmfld(1)*zdif

            fzz(i) = fzz(i) + gamma
            engfld = engfld - gamma*zdif/2.

          endif

        enddo

      else

c     unidentified field potential error exit

        call error(idnode,454)

      endif

#ifdef VAMPIR
      call VTEND(27, ierr)
#endif
      return
      end
