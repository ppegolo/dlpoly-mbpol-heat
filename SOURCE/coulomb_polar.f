      subroutine coulomb_polar
     x    (idnode,imcon,mxnode,natms,engcpe,vircpe,
     x    efieldkx,efieldky,efieldkz,
     x    efdcrecx,efdcrecy,efdcrecz,
     x    efddmurecx,efddmurecy,efddmurecz,
     x    dipx,dipy,dipz,emux,emuy,emuz,polr2)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating coulombic energy in a
c     polar system from charge-dipole and dipole-dipole interactions
c
c     parallel replicated data version
c
c     copyright - voth group
c     author    - c. j. burnham
c
c     2003/10/04
c
c***********************************************************************
c
#ifdef HEAT_CURRENT
      use heatcurrent, only: update_energy_polar
#endif

#include "dl_params.inc"

      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension efdcrecx(mxatms),efdcrecy(mxatms),efdcrecz(mxatms)
      dimension efddmurecx(mxatms),efddmurecy(mxatms),
     x          efddmurecz(mxatms)
      dimension polr2(mxatms)


#ifdef VAMPIR
      call VTBEGIN(100, ierr)
#endif
c
c     no external field

      extx=0.d0
      exty=0.d0
      extz=0.d0

      engcpe=0.d0
      vircpe=0.d0

      ucd=0.d0
      udd=0.d0
      uspring=0.d0

c     set up atoms numbers for nodes
      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode

      do i=iatm0,iatm1

         if (polr2(i).gt.1.d-6) then

            engcpe = engcpe-
     x           (dipx(i)*(efieldkx(i)+efdcrecx(i)+
     x            0.5d0*(emux(i)+efddmurecx(i))+extx)+
     x            dipy(i)*(efieldky(i)+efdcrecy(i)+
     x            0.5d0*(emuy(i)+efddmurecy(i))+exty)+
     x            dipz(i)*(efieldkz(i)+efdcrecz(i)+
     x            0.5d0*(emuz(i)+efddmurecz(i))+extz))

            ucd=ucd-
     x           (dipx(i)*(efieldkx(i)+efdcrecx(i))+
     x            dipy(i)*(efieldky(i)+efdcrecy(i))+
     x            dipz(i)*(efieldkz(i)+efdcrecz(i)))

            udd=udd-
     x        (dipx(i)*0.5d0*(emux(i)+efddmurecx(i))+
     x         dipy(i)*0.5d0*(emuy(i)+efddmurecy(i))+
     x         dipz(i)*0.5d0*(emuz(i)+efddmurecz(i)))
c
c     energy to induce dipole

            edum = dipx(i)**2+dipy(i)**2+dipz(i)**2
            engcpe = engcpe + r4pie0*edum/polr2(i)/2.d0
            uspring=uspring+r4pie0*edum/polr2(i)/2.d0

#ifdef HEAT_CURRENT
            call update_energy_polar(i,
     x      r4pie0*edum/polr2(i)/2.d0)
#endif

         endif

      enddo

#if 0
c        engcpe=0.d0
      if (mxnode.gt.1) then
         call gdsum(ucd,1,edum)
         call gdsum(udd,1,edum)
         call gdsum(uspring,1,edum)
      end if

      if (idnode.eq.0) then
        write(50,*)'uspring =',uspring/418.4d0
        write(50,*)'ucd =',0.5*ucd/418.4d0
        write(50,*)'udd =',udd/418.4d0
      end if

      scale = 18.22234397655801030455d0 ! VB: for comparison with c++
      do i = iatm0, iatm1
          write(50+idnode,'(i2,3(1x,ES12.5))') i,
     x scale*dipx(i), scale*dipy(i), scale*dipz(i)
      end do
#endif

#ifdef VAMPIR
      call VTEND(100, ierr)
#endif
      return
      end
