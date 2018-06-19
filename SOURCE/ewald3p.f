! 30 NOV 05 - IUCHI - ESP FOR DMS
! 10 NOV 05 - IUCHI - INTRODUCING FORCE DECOMPOSITION
! 10 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
!
      subroutine ewald3p
     x   (iatm,ilst,engcpe,vircpe,alpha,epsq,nexatm,
     x   lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x   dipx,dipy,dipz,polr,polr2,athole12,athole13,nthole,lthole,
     x   lttm,nttm2,listttm2,n_water,natms)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method. This is the subrotine
c     for correction of the excluded pairs.
c
c     parallel replicated data version (part 3)
c
c     copyright - daresbury laboratory 1992
c                 voth group
c     author    - w. smith dec 1992.
c                 c. j. burnham  sept 2003.
c                 t. yan sept 2003.
c
c     stress stensor added t.forester may 1994
!
!     Last updated: 30 Nov 2005 by S. Iuchi
c
c***********************************************************************
c
! from MODULE
#ifdef TTM_FORCE_DECOMPOSITION
      use ttm_forces,  only: fx_ttm_per, fy_ttm_per, fz_ttm_per,
     $                       fx_ttm_ind, fy_ttm_ind, fz_ttm_ind
#endif /* TTM_FORCE_DECOMPOSITION */

#ifdef HEAT_CURRENT
      use heatcurrent, only: update_stress_ew3p, update_energy_ew3p,
     x                       update_forces
#endif /* HEAT_CURRENT */

      use ps_type_dms, only: vesp_c, vesp_dc_c, ldms

#include "dl_params.inc"

      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension chge(mxatms),polr(mxatms),polr2(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension listttm2(mxatms)
      dimension stress(9)

      dimension bbb(0:3),chgchg(0:3),chgdip(0:3),dipdip(0:3)

      logical lthole,lttm

      integer n_water, n_water_atoms, no_water_atoms
#ifdef HEAT_CURRENT
#ifdef HEAT_STRESS
      real(8) :: strs1,strs2,strs3,strs4,strs5,strs6
      real(8) :: strs7,strs8,strs9
#endif
      real(8) :: force_tmp(3)
#endif /* HEAT_CURRENT */

#ifdef VAMPIR
      call VTBEGIN(103, ierr)
#endif

c
c     initialise potential energy and virial

      engcpe=0.d0
      vircpe=0.d0
      espdum=0.d0    ! for ESP
      espdcdum=0.0d0 ! for ESP

c
c     start of primary loop for forces evaluation

      do m=1,nexatm(ilst)

c
c     atomic index and charge product

        jatm=lexatm(ilst,m)
c
c     calculate interatomic distance

        rsq=xdf(m)**2+ydf(m)**2+zdf(m)**2
        rrr = dsqrt(rsq)

        rri = 1.d0/rrr
        rsqi = 1.d0/rsq
        r3i = 1.d0/(rrr**3)
        r5i = r3i*rsqi
c
c     calculate the 'Smith' B (error) functions and gradients
c     for real space Ewald sum.  (I've modified
c     Smith's recursion to give -erf, not erfc functions.)
c     comment by c. j. burnham


        a=alpha*rrr
        bbb(0)=(erfcc(a)-1.d0)*rri
        exp2a=exp(-a*a)

        do n=1,3

           fn=dble(n)
           bbb(n)=rsqi*((2.0*fn-1.d0)*bbb(n-1)
     x         +((2.0*alpha**2.0)**fn)*exp2a/alpha/sqrpi)

        enddo
c
c     assume same screening functions for chgs and dips

        do n=0,3

          chgchg(n)=bbb(n)
          chgdip(n)=bbb(n)
          dipdip(n)=bbb(n)

        enddo
c
c     if ewald parameter is too small then turn it off!

              if(alpha.lt.1.d-6) then
                do n=0,3
                  chgchg(n)=0.d0
                  chgdip(n)=0.d0
                  dipdip(n)=0.d0
                enddo
              endif
c
c     remove charge-charge interactions

           cicj=chge(iatm)*chge(jatm)

           engcpe = engcpe+cicj*chgchg(0)*r4pie0/epsq
#ifdef HEAT_CURRENT
      call update_energy_ew3p(iatm,0.5d0*cicj*chgchg(0)*
     x       r4pie0/epsq)
      call update_energy_ew3p(jatm,0.5d0*cicj*chgchg(0)*
     x       r4pie0/epsq)
#endif /* HEAT_CURRENT */
!
!     for ESP
           espdum   = chgchg(0) * r4pie0 / epsq
           espdcdum = chgdip(1) * r4pie0 / epsq

           ffac = cicj*chgchg(1)*r4pie0/epsq

           ffac_cc = cicj*chgchg(1)*r4pie0/epsq  ! for force decompositions
c
c     remove charge-dipole interactions

           rdmui = dipx(iatm)*xdf(m)+dipy(iatm)*ydf(m)+
     x             dipz(iatm)*zdf(m)
           rdmuj = dipx(jatm)*xdf(m)+dipy(jatm)*ydf(m)+
     x             dipz(jatm)*zdf(m)

           facmu = (chge(jatm)*rdmui-chge(iatm)*rdmuj)*
     x             (-chgdip(2))*r4pie0/epsq

           dfacmui =  chgdip(1)*chge(jatm)*r4pie0/epsq
           dfacmuj = -chgdip(1)*chge(iatm)*r4pie0/epsq
c
c     add intramolecular dipole-dipole interactions
c
           if (polr2(iatm).gt.1.d-6 .and. polr2(jatm).gt.1.d-6) then

             dpidpj = dipx(iatm)*dipx(jatm)+dipy(iatm)*dipy(jatm)+
     x                dipz(iatm)*dipz(jatm)

             facmu = facmu + (dipdip(2)*dpidpj-
     x               dipdip(3)*rdmui*rdmuj)*r4pie0/epsq

             dfacmui = dfacmui+dipdip(2)*rdmuj*r4pie0/epsq
             dfacmuj = dfacmuj+dipdip(2)*rdmui*r4pie0/epsq
c
c    add Applequist's point dipole model

             facmu1 = -15.d0*rdmui*rdmuj*rsqi*r5i*r4pie0/epsq
             facmu2 =   3.d0*dpidpj*r5i*r4pie0/epsq
             facmu  = facmu+facmu1+facmu2

             dfacmuidd = 3.d0*rdmuj*r5i*r4pie0/epsq
             dfacmujdd = 3.d0*rdmui*r5i*r4pie0/epsq

             dfacmui = dfacmui + dfacmuidd
             dfacmuj = dfacmuj + dfacmujdd
c
c     Thole's attenuated dipole model

             if (lthole) then

               bb=(polr(iatm)*polr(jatm))**(1.d0/6.d0)

!FP_fix_start
               n_water_atoms = 4 * n_water
               no_water_atoms = natms - n_water_atoms

               athole=athole13
               if (mod((iatm-no_water_atoms)-1,3) .eq. 0
     x          .or. mod((jatm-no_water_atoms)-1,3) .eq. 0)
     x             athole=athole12

c              write(910,'(4i8)') natms, no_water_atoms,
c    x                            n_water, n_water_atoms
c              write(911,'(2i8,5f12.5)') iatm, jatm,
c    x            polr2(iatm), polr2(jatm),
c    x            athole, athole12, athole13
c
c              athole=athole13
c              if(mod(iatm-1,3).eq.0.or.mod(jatm-1,3).eq.0)
c    x             athole=athole12
c              write(910,'(4i8)') natms, no_water_atoms,
c    x                            n_water, n_water_atoms
c              write(911,'(2i8,5f12.5)') iatm, jatm,
c    x            polr2(iatm), polr2(jatm),
c    x            athole, athole12, athole13
!FP_fix_end

               if (nthole.eq.3) then

                 arob3=athole*(rrr/bb)**3
                 dum1=dexp(-arob3)
                 dum2=(1.d0+arob3)*dum1
                 dum3=9.d0*athole**2*dum1/bb**6
                 dum4=3.d0*athole*dum1/bb**3

                 facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                   rdmui*rdmuj*dum3*rri*r4pie0/epsq-
     x                   dpidpj*rri*rri*dum4*r4pie0/epsq

               else if (nthole.eq.4) then

                 arob4=athole*(rrr/bb)**4
                 dum1=dexp(-arob4)
                 dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
                 dum3=4.d0*athole*(4.d0*arob4-1.d0)*dum1*rsqi/bb**4
                 dum4=4.d0*athole*dum1*rrr/bb**4

                 facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                   rdmui*rdmuj*dum3*rri*r4pie0/epsq-
     x                   dpidpj*rsqi*dum4*r4pie0/epsq

               else

                 arob  = athole*rrr/bb
                 arob2 = arob*arob
                 arob3 = arob2*arob

                 dum1 = dexp(-arob)*(arob2/2.d0+arob+1.d0)
                 dum2 = dum1+dexp(-arob)*arob3/6.d0
                 dum3 = (athole/bb)*dexp(-arob)/2.d0
                 dum4 = dum3*arob2
                 dum5 = dum3*arob3

                 facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                   rdmui*rdmuj*rri*r5i*dum5*r4pie0/epsq-
     x                   dpidpj*rri*r3i*dum4*r4pie0/epsq

               endif

               dfacmui = dfacmui - dfacmuidd*dum2
               dfacmuj = dfacmuj - dfacmujdd*dum2

             endif

           endif ! polr2(iatm).gt.1.d-6 .and. polr2(jatm).gt.1.d-6

           dfacx = dfacmui*dipx(iatm)+dfacmuj*dipx(jatm)
           dfacy = dfacmui*dipy(iatm)+dfacmuj*dipy(jatm)
           dfacz = dfacmui*dipz(iatm)+dfacmuj*dipz(jatm)
c
c     increment forces

           ffac = ffac + facmu

           fx = ffac*xdf(m)+dfacx
           fy = ffac*ydf(m)+dfacy
           fz = ffac*zdf(m)+dfacz

           fxx(iatm) = fxx(iatm)+fx
           fyy(iatm) = fyy(iatm)+fy
           fzz(iatm) = fzz(iatm)+fz

           fxx(jatm) = fxx(jatm)-fx
           fyy(jatm) = fyy(jatm)-fy
           fzz(jatm) = fzz(jatm)-fz

#ifdef HEAT_CURRENT
            force_tmp = (/ fx,fy,fz /)
            call update_forces(iatm,jatm,force_tmp)
            call update_forces(jatm,iatm,-force_tmp)
#endif /*HEAT_CURRENT*/

           virdum = dfacx*xdf(m)+dfacy*ydf(m)+dfacz*zdf(m)
           vircpe = vircpe - ffac*rsq - virdum
!
!     for force decompositions
#ifdef TTM_FORCE_DECOMPOSITION
           fx_cc = ffac_cc * xdf(m)
           fy_cc = ffac_cc * ydf(m)
           fz_cc = ffac_cc * zdf(m)

           fx_dc_dd = facmu * xdf(m) + dfacx
           fy_dc_dd = facmu * ydf(m) + dfacy
           fz_dc_dd = facmu * zdf(m) + dfacz

           fx_ttm_per(iatm) = fx_ttm_per(iatm) + fx_cc
           fy_ttm_per(iatm) = fy_ttm_per(iatm) + fy_cc
           fz_ttm_per(iatm) = fz_ttm_per(iatm) + fz_cc

           fx_ttm_per(jatm) = fx_ttm_per(jatm) - fx_cc
           fy_ttm_per(jatm) = fy_ttm_per(jatm) - fy_cc
           fz_ttm_per(jatm) = fz_ttm_per(jatm) - fz_cc

           fx_ttm_ind(iatm) = fx_ttm_ind(iatm) + fx_dc_dd
           fy_ttm_ind(iatm) = fy_ttm_ind(iatm) + fy_dc_dd
           fz_ttm_ind(iatm) = fz_ttm_ind(iatm) + fz_dc_dd

           fx_ttm_ind(jatm) = fx_ttm_ind(jatm) - fx_dc_dd
           fy_ttm_ind(jatm) = fy_ttm_ind(jatm) - fy_dc_dd
           fz_ttm_ind(jatm) = fz_ttm_ind(jatm) - fz_dc_dd
#endif /* TTM_FORCE_DECOMPOSITION */

           if( ldms ) then  ! ESP for DMS

              vesp_c(iatm) = vesp_c(iatm) + chge(jatm) * espdum
              vesp_c(jatm) = vesp_c(jatm) + chge(iatm) * espdum

              vesp_c(iatm) = vesp_c(iatm) + rdmuj * espdcdum
              vesp_c(jatm) = vesp_c(jatm) - rdmui * espdcdum

           end if

#ifdef STRESS
c
c     calculate stress tensor

        stress(1) = stress(1)+ xdf(m)*fx
        stress(2) = stress(2)+ xdf(m)*fy
        stress(3) = stress(3)+ xdf(m)*fz

        stress(5) = stress(5)+ ydf(m)*fy
        stress(6) = stress(6)+ ydf(m)*fz

        stress(9) = stress(9)+ zdf(m)*fz
#endif

#ifdef HEAT_CURRENT
#ifdef HEAT_STRESS
        strs1 = strs1+xdf(m)*fx
        strs2 = strs2+xdf(m)*fy
        strs3 = strs3+xdf(m)*fz
        strs4 = strs4+xdf(m)*fy
        strs5 = strs5+ydf(m)*fy
        strs6 = strs6+ydf(m)*fz
        strs7 = strs7+xdf(m)*fz
        strs8 = strs8+ydf(m)*fz
        strs9 = strs9+zdf(m)*fz

        call update_stress_ew3p(iatm,1,1,-0.5d0*xdf(m)*fx)
        call update_stress_ew3p(iatm,1,2,-0.5d0*xdf(m)*fy)
        call update_stress_ew3p(iatm,1,3,-0.5d0*xdf(m)*fz)
        call update_stress_ew3p(iatm,2,1,-0.5d0*xdf(m)*fy)
        call update_stress_ew3p(iatm,2,2,-0.5d0*ydf(m)*fy)
        call update_stress_ew3p(iatm,2,3,-0.5d0*ydf(m)*fz)
        call update_stress_ew3p(iatm,3,1,-0.5d0*xdf(m)*fz)
        call update_stress_ew3p(iatm,3,2,-0.5d0*ydf(m)*fz)
        call update_stress_ew3p(iatm,3,3,-0.5d0*zdf(m)*fz)
        call update_stress_ew3p(jatm,1,1,-0.5d0*xdf(m)*fx)
        call update_stress_ew3p(jatm,1,2,-0.5d0*xdf(m)*fy)
        call update_stress_ew3p(jatm,1,3,-0.5d0*xdf(m)*fz)
        call update_stress_ew3p(jatm,2,1,-0.5d0*xdf(m)*fy)
        call update_stress_ew3p(jatm,2,2,-0.5d0*ydf(m)*fy)
        call update_stress_ew3p(jatm,2,3,-0.5d0*ydf(m)*fz)
        call update_stress_ew3p(jatm,3,1,-0.5d0*xdf(m)*fz)
        call update_stress_ew3p(jatm,3,2,-0.5d0*ydf(m)*fz)
        call update_stress_ew3p(jatm,3,3,-0.5d0*zdf(m)*fz)
#endif
#endif /* HEAT_CURRENT */
      enddo
#ifdef STRESS
c
c     complete stress tensor

      stress(4) = stress(2)
      stress(7) = stress(3)
      stress(8) = stress(6)
#endif

#ifdef VAMPIR
      call VTEND(103, ierr)
#endif
      return
      end
