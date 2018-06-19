#include "assert.h"
! 30 NOV 05 - IUCHI - ESP FOR DMS
! 10 NOV 05 - IUCHI - INTRODUCING FORCE DECOMPOSITION
! 10 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
!
      subroutine ewald2p
     x  (iatm,ik,engcpe,vircpe,rcut,epsq,ilist,
     x  chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x  polr,polr2,athole,athole_ion,ithole,n_ions,
     x  athole_ionwat,
     x  nthole,lthole,lacs,lads,
     x  ascc,ascd,dipx,dipy,dipz,alpha,ercp,drewd,
     x  lttm,nttm2,listttm2)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating electric energy in a
c     periodic system using ewald's method
c
c     parallel replicated data version (part 2)
c
c     copyright - voth group
c     author    - c. j. burnham
c     parallel  - t. yan
c
c     part 2 - real space terms.
c
c     Tabulated potential in r space
c
c     2003/10/01 19:47:27
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
      use ps_type_dms, only: vesp_r, vesp_dc_r, ldms

#ifdef HEAT_CURRENT
      use heatcurrent, only: update_stress_ew2p, update_energy_ew2p,
     x                       update_forces
#endif

#include "dl_params.inc"

      parameter(g23=1.3541179394264d0)

      logical lthole,lacs,lads,lttm

      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ilist(mxxdf)
      dimension chge(mxatms),polr(mxatms),polr2(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension listttm2(mxatms)
      dimension ercp(mxegrd,0:3)
      dimension stress(9)
      dimension fi(3)
#ifdef TTM_FORCE_DECOMPOSITION
      dimension fi_cc(3), fi_dc_dd(3)  ! for force decompositions
#endif /* TTM_FORCE_DECOMPOSITION  */

      dimension bbb(0:3),chgchg(0:3),chgdip(0:3),dipdip(0:3)

!FP_fix_start
      real(8) athole_ion,athole_ionwat
      integer n_water, n_water_atoms, no_water_atoms, n_ions, ithole
!FP_fix_end

#ifdef HEAT_CURRENT
      real(8), parameter :: third=1.0d0/3.0d0
      real(8) :: force_tmp(3)
#endif


#ifdef VAMPIR
      call VTBEGIN(100, ierr)
#endif
c
c     set cutoff condition for pair forces

      rcsq=rcut**2
c
c     reciprocal of interpolation interval

#ifdef STRESS
c
c     initialise stress tensor accumulators
      strs1 = 0.d0
      strs2 = 0.d0
      strs3 = 0.d0
      strs5 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0
#endif
c
c     external field is zero

      extx=0.d0
      exty=0.d0
      extz=0.d0
c
c     initialise potential energy and virial

      engcpe=0.d0
      vircpe=0.d0

c
c     start of primary loop for forces evaluation

ccc      chgea = chge(iatm)/epsq*r4pie0
ccc      if(chgea.ne.0.d0) then

        fi(1) = fxx(iatm)
        fi(2) = fyy(iatm)
        fi(3) = fzz(iatm)

!     for force decompositions
#ifdef TTM_FORCE_DECOMPOSITION
        fi_cc(1) = fx_ttm_per(iatm)
        fi_cc(2) = fy_ttm_per(iatm)
        fi_cc(3) = fz_ttm_per(iatm)

        fi_dc_dd(1) = fx_ttm_ind(iatm)
        fi_dc_dd(2) = fy_ttm_ind(iatm)
        fi_dc_dd(3) = fz_ttm_ind(iatm)
#endif /* TTM_FORCE_DECOMPOSITION */

        do m=1,ik
c
c     atomic index and charge product

          jatm=ilist(m)
!
! for ESP

          espdum=0.0d0
          espdcdum=0.0d0
c
c     calculate interatomic distance

            rsq=rsqdf(m)
c
c     apply truncation of potential

            if(rcsq.gt.rsq)then

c
c     calculate interatomic distance

              rrr = sqrt(rsq)

              ll = int(rrr/drewd)
              l1=ll+1
              l2=ll+2
              ppp = rrr/drewd - dble(ll)
c
c     calculate interaction energy using 3-point interpolation

              do n=0,3

                 vk = ercp(ll,n)
                 vk1 = ercp(l1,n)
                 vk2 = ercp(l2,n)

                 t1 = vk + (vk1-vk)*ppp
                 t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)

                 bbb(n) = t1 + (t2-t1)*ppp*0.5d0

              enddo

              rri = 1.d0/rrr
              rsqi = 1.d0/rsq
              r3i = 1.d0/(rrr**3)
              r5i = r3i*rsqi
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
c     no screening or switching functions for charge and dipole

              screencc=1.d0
              dscreencc=0.d0
              screendc=1.d0
              dscreendc=0.d0
              screendd=1.d0
              dscreendd=0.d0
c
c     charge-charge interactions

              cicj=chge(iatm)*chge(jatm)

              edum = cicj*(chgchg(0)+screencc/rrr)*r4pie0
              engcpe=engcpe+edum

#ifdef HEAT_CURRENT
              call update_energy_ew2p(iatm,0.50d0*edum)
              call update_energy_ew2p(jatm,0.5d0*edum)
#endif /* HEAT_CURRENT */

              espdum = (chgchg(0)+screencc/rrr)*r4pie0  ! ESP

              ffac = cicj*(screencc*r3i-dscreencc*rsqi+
     x               chgchg(1))*r4pie0

!     for force decompositions
#ifdef TTM_FORCE_DECOMPOSITION
              ffac_cc = cicj*(screencc*r3i-dscreencc*rsqi+
     x             chgchg(1))*r4pie0
#endif /* TTM_FORCE_DECOMPOSITION */


              if (lthole .and. lacs) then
c
c     Smeared charge-charge interactions

                bb=(polr(iatm)*polr(jatm))**(1.d0/6.d0)

                if (nthole.eq.3) then
                  arob3=ascc*(rrr/bb)**3
                  dum1=dexp(-arob3)
                  engcpe=engcpe-cicj*(dum1-arob3**(1.d0/3.d0)
     x                  *g23*gammq(2.d0/3.d0,arob3))/rrr*r4pie0

#ifdef HEAT_CURRENT
       call update_energy_ew2p(iatm,-cicj*0.5d0*(dum1-arob3**third*
     x                  g23*gammq(2.d0/3.d0,arob3))/rrr*r4pie0)
       call update_energy_ew2p(jatm,-cicj*0.5d0*(dum1-arob3**third*
     x                  g23*gammq(2.d0/3.d0,arob3))/rrr*r4pie0)
#endif /* HEAT_CURRENT */


                  ffac=ffac-cicj*r3i*dum1*r4pie0
!     for ESP
                  espdum=espdum-(dum1-arob3**(1.d0/3.d0)
     x                  *g23*gammq(2.d0/3.d0,arob3))/rrr*r4pie0
                else if (nthole.eq.4) then
                  arob4=ascc*(rrr/bb)**4
                  dum1=dexp(-arob4)
                  g34=dexp(gammln(3.d0/4.d0))
                  engcpe=engcpe-cicj*(dum1-arob4**(1.d0/4.d0)
     x                  *g34*gammq(3.d0/4.d0,arob4))/rrr*r4pie0

#ifdef HEAT_CURRENT
       call update_energy_ew2p(iatm,-cicj*0.5d0*(dum1-arob4**0.25d0
     x                  *g34*gammq(0.75d0,arob4))/rrr*r4pie0)
       call update_energy_ew2p(jatm,-cicj*0.5d0*(dum1-arob4**0.25d0
     x                  *g34*gammq(0.75d0,arob4))/rrr*r4pie0)
#endif /* HEAT_CURRENT */

                  espdum=espdum-(dum1-arob4**(1.d0/4.d0)
     x                  *g34*gammq(3.d0/4.d0,arob4))/rrr*r4pie0
                  ffac=ffac-cicj*r3i*dum1*r4pie0
                else
                  arob=ascc*rrr/bb
                  arob2=arob*arob
                  dum1=dexp(-arob)*(arob2/2.d0+arob+1.d0)
                  engcpe=engcpe-cicj*dum1/rrr*r4pie0+
     x               0.5d0*cicj*(ascc/bb)*dexp(-arob)*(1+arob)*r4pie0

#ifdef HEAT_CURRENT
       call update_energy_ew2p(iatm,0.5d0*(-cicj*dum1/rrr*r4pie0+
     x              0.5d0*cicj*(ascc/bb)*dexp(-arob)*(1+arob)*r4pie0))
       call update_energy_ew2p(jatm,0.5d0*(-cicj*dum1/rrr*r4pie0+
     x              0.5d0*cicj*(ascc/bb)*dexp(-arob)*(1+arob)*r4pie0))
#endif /* HEAT_CURRENT */

                  ffac=ffac-cicj*r3i*dum1*r4pie0
                endif

!     for force decompositions
#ifdef TTM_FORCE_DECOMPOSITION
                   ffac_cc=ffac_cc-cicj*r3i*dum1*r4pie0
#endif /* TTM_FORCE_DECOMPOSITION */

              endif ! lthole .and. lacs

c
c     charge-dipole

              rdmui = dipx(iatm)*xdf(m)+dipy(iatm)*ydf(m)+
     x                dipz(iatm)*zdf(m)
              rdmuj = dipx(jatm)*xdf(m)+dipy(jatm)*ydf(m)+
     x                dipz(jatm)*zdf(m)

              facmu0 = (chge(jatm)*rdmui-chge(iatm)*rdmuj)*
     x                 (rrr*dscreendc-3.d0*screendc)*r5i*r4pie0
              facmui0 = screendc*r3i*chge(jatm)*r4pie0
              facmuj0 =-screendc*r3i*chge(iatm)*r4pie0

              facmu = facmu0 - (chge(jatm)*rdmui-chge(iatm)*rdmuj)*
     x                          chgdip(2)*r4pie0

              dfacmui = facmui0+chgdip(1)*chge(jatm)*r4pie0
              dfacmuj = facmuj0-chgdip(1)*chge(iatm)*r4pie0
!
!     for ESP (charge-dipole)

              espdcdum = ( chgdip(1) + r3i ) * r4pie0
c
c     option for the all dipole smeared, i.e. lads=.true.

              if (lthole .and. lads) then

                bb = (polr(iatm)*polr(jatm))**(1.d0/6.d0)

                assert(bb.gt.1.d-6)

                if (nthole.eq.3) then

                  arob3=ascd*(rrr/bb)**3
                  dum1=dexp(-arob3)
                  dum2=(1.d0+arob3)*dum1

                else if (nthole.eq.4) then

                  arob4=ascd*(rrr/bb)**4
                  dum1=dexp(-arob4)
                  dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1

                else

                  arob  = ascd*rrr/bb
                  arob2 = arob*arob
                  arob3 = arob2*arob

                  dum1 = dexp(-arob)*(arob2/2.d0+arob+1.d0)
                  dum2 = dum1+dexp(-arob)*arob3/6.d0

                endif
c
c     smearing charge-dipole

                facmu=facmu-facmu0*dum2

c
c     smearing charge-dipole wrt dipole direction

                dfacmui = dfacmui - facmui0*dum1
                dfacmuj = dfacmuj - facmuj0*dum1
!
!     for ESP (charge-dipole)

                espdcdum = espdcdum - r3i * dum1 * r4pie0

              endif
c
c     Ignore interaction if atoms' polarizabilites are zero

          if (polr2(iatm).gt.1.d-6 .and. polr2(jatm).gt.1.d-6) then
c
c     dipole-dipole interactions
!             write(5002,'(2i8,2f15.5)') iatm,jatm,dipx(iatm),dipx(jatm)

              dpidpj = dipx(iatm)*dipx(jatm)+dipy(iatm)*dipy(jatm)+
     x                 dipz(iatm)*dipz(jatm)

              facmu = facmu+((dpidpj*(3.d0*screendd-dscreendd/rrr)+
     x                rdmui*rdmuj*rsqi*(3.d0*dscreendd*rrr-
     x                15.d0*screendd))*r5i+
     x                dipdip(2)*dpidpj-dipdip(3)*rdmui*rdmuj)*r4pie0

              dfacmui=dfacmui+(screendd*3.d0*r5i+dipdip(2))*rdmuj*r4pie0
              dfacmuj=dfacmuj+(screendd*3.d0*r5i+dipdip(2))*rdmui*r4pie0
c
c     option for the all dipole smeared, i.e. lads=.true.

              if (lthole .and. lads) then

                facmu1 = -15.d0*rdmui*rdmuj*rsqi*r5i*r4pie0
                facmu2 =   3.d0*dpidpj*r5i*r4pie0

                dfacmui1 = 3.d0*rdmuj*r5i*r4pie0
                dfacmuj1 = 3.d0*rdmui*r5i*r4pie0

                bb = (polr(iatm)*polr(jatm))**(1.d0/6.d0)

                if (nthole.eq.3) then

                  arob3=athole*(rrr/bb)**3
                  dum1=dexp(-arob3)
                  dum2=(1.d0+arob3)*dum1
                  dum3=9.d0*athole**2*dum1/bb**6
                  dum4=3.d0*athole*dum1/bb**3

                  facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                    rdmui*rdmuj*dum3*rri*r4pie0-
     x                    dpidpj*rri*rri*dum4*r4pie0

                else if (nthole.eq.4) then
!FP_fix_start: different Thole damping and smearing for ion-ion.
                  if ((iatm.le.n_ions) .and. (jatm.le.n_ions)) then
                    if (ithole.eq.3) then
                      arob3=athole_ion*(rrr/bb)**3
                      dum1=dexp(-arob3)
                      dum2=(1.d0+arob3)*dum1
                      dum3=9.d0*athole_ion**2*dum1/bb**6
                      dum4=3.d0*athole_ion*dum1/bb**3

                      facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                        rdmui*rdmuj*dum3*rri*r4pie0-
     x                        dpidpj*rri*rri*dum4*r4pie0
c                     write(901,'(4i8,5f12.5)')
c    x                       iatm, jatm, n_ions, ithole,
c    x                       polr(iatm), polr(jatm), athole_ion
                    else
                      arob4=athole_ion*(rrr/bb)**4
                      dum1=dexp(-arob4)
                      dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
                      dum3=4.d0*athole_ion
     x                    *(4.d0*arob4-1.d0)*dum1*rsqi/bb**4
                      dum4=4.d0*athole_ion*dum1*rrr/bb**4

                      facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                        rdmui*rdmuj*dum3*rri*r4pie0-
     x                        dpidpj*rsqi*dum4*r4pie0
c                     write(902,'(4i8,5f12.5)')
c    x                       iatm, jatm, n_ions, ithole,
c    x                       polr(iatm), polr(jatm), athole_ion
                    end if
c  for ion-water

        elseif (((iatm.le.n_ions).and.(jatm.gt.n_ions))
     x        .or.((iatm.gt.n_ions).and.(jatm.le.n_ions))) then

                      arob4=athole_ionwat*(rrr/bb)**4
                      dum1=dexp(-arob4)
                      dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
                      dum3=4.d0*athole_ionwat
     x                    *(4.d0*arob4-1.d0)*dum1*rsqi/bb**4
                      dum4=4.d0*athole_ionwat*dum1*rrr/bb**4

                      facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                        rdmui*rdmuj*dum3*rri*r4pie0-
     x                        dpidpj*rsqi*dum4*r4pie0

c                      write(4252,*) iatm,jatm

                  else
                    arob4=athole*(rrr/bb)**4
                    dum1=dexp(-arob4)
                    dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
                    dum3=4.d0*athole*(4.d0*arob4-1.d0)*dum1*rsqi/bb**4
                    dum4=4.d0*athole*dum1*rrr/bb**4

                    facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                      rdmui*rdmuj*dum3*rri*r4pie0-
     x                      dpidpj*rsqi*dum4*r4pie0
c                      write(4252,*) iatm,jatm
c                   write(903,'(2i8,5f12.5)')
c    x                    iatm, jatm, polr(iatm), polr(jatm), athole
                  end if
!FP_fix_end: different Thole damping and smearing for ion-ion.
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
     x                    rdmui*rdmuj*rri*r5i*dum5*r4pie0-
     x                    dpidpj*rri*r3i*dum4*r4pie0
                endif
c
c     smearing dipole-dipole wrt dipole direction

                dfacmui = dfacmui - dfacmui1*dum2
                dfacmuj = dfacmuj - dfacmuj1*dum2

              endif

              endif ! polr2(iatm).gt.1.d-6 .and. polr2(jatm).gt.1.d-6

              dfacx = dfacmui*dipx(iatm)+dfacmuj*dipx(jatm)
              dfacy = dfacmui*dipy(iatm)+dfacmuj*dipy(jatm)
              dfacz = dfacmui*dipz(iatm)+dfacmuj*dipz(jatm)
c
c     increment forces

              ffac = ffac+facmu

              virdum = dfacx*xdf(m)+dfacy*ydf(m)+dfacz*zdf(m)
              vircpe = vircpe - ffac*rsq - virdum

              fx = ffac*xdf(m)+dfacx
              fy = ffac*ydf(m)+dfacy
              fz = ffac*zdf(m)+dfacz

              fi(1) = fi(1)+fx
              fi(2) = fi(2)+fy
              fi(3) = fi(3)+fz

              fxx(jatm) = fxx(jatm)-fx
              fyy(jatm) = fyy(jatm)-fy
              fzz(jatm) = fzz(jatm)-fz

#ifdef HEAT_CURRENT
              force_tmp = (/-fx,-fy,-fz/)
              call update_forces(jatm,iatm,force_tmp)
              call update_forces(iatm,jatm,-force_tmp)
#endif /*HEAT_CURRENT*/

!     for force decompositions
#ifdef TTM_FORCE_DECOMPOSITION
              fx_cc = ffac_cc * xdf(m)
              fy_cc = ffac_cc * ydf(m)
              fz_cc = ffac_cc * zdf(m)

              fi_cc(1) = fi_cc(1) + fx_cc
              fi_cc(2) = fi_cc(2) + fy_cc
              fi_cc(3) = fi_cc(3) + fz_cc

              fx_ttm_per(jatm) = fx_ttm_per(jatm) - fx_cc
              fy_ttm_per(jatm) = fy_ttm_per(jatm) - fy_cc
              fz_ttm_per(jatm) = fz_ttm_per(jatm) - fz_cc

              fx_dc_dd = facmu * xdf(m) + dfacx
              fy_dc_dd = facmu * ydf(m) + dfacy
              fz_dc_dd = facmu * zdf(m) + dfacz

              fi_dc_dd(1) = fi_dc_dd(1) + fx_dc_dd
              fi_dc_dd(2) = fi_dc_dd(2) + fy_dc_dd
              fi_dc_dd(3) = fi_dc_dd(3) + fz_dc_dd

              fx_ttm_ind(jatm) = fx_ttm_ind(jatm) - fx_dc_dd
              fy_ttm_ind(jatm) = fy_ttm_ind(jatm) - fy_dc_dd
              fz_ttm_ind(jatm) = fz_ttm_ind(jatm) - fz_dc_dd
#endif /* TTM_FORCE_DECOMPOSITION */
!     ESP for DMS

              if( ldms ) then

                 vesp_r(iatm) = vesp_r(iatm) + chge(jatm) * espdum
                 vesp_r(jatm) = vesp_r(jatm) + chge(iatm) * espdum

                 vesp_dc_r(iatm) = vesp_dc_r(iatm) + rdmuj * espdcdum
                 vesp_dc_r(jatm) = vesp_dc_r(jatm) - rdmui * espdcdum

              endif

#ifdef STRESS
c
c     calculate stress tensor

              strs1 = strs1 + xdf(m)*fx
              strs2 = strs2 + xdf(m)*fy
              strs3 = strs3 + xdf(m)*fz

              strs5 = strs5 + ydf(m)*fy
              strs6 = strs6 + ydf(m)*fz

              strs9 = strs9 + zdf(m)*fz

#endif

#ifdef HEAT_CURRENT
#ifdef HEAT_STRESS
        call update_stress_ew2p(jatm,1,1,-0.5d0*xdf(m)*fx)
        call update_stress_ew2p(jatm,1,2,-0.5d0*xdf(m)*fy)
        call update_stress_ew2p(jatm,1,3,-0.5d0*xdf(m)*fz)
        call update_stress_ew2p(jatm,2,1,-0.5d0*xdf(m)*fy)
        call update_stress_ew2p(jatm,2,2,-0.5d0*ydf(m)*fy)
        call update_stress_ew2p(jatm,2,3,-0.5d0*ydf(m)*fz)
        call update_stress_ew2p(jatm,3,1,-0.5d0*xdf(m)*fz)
        call update_stress_ew2p(jatm,3,2,-0.5d0*ydf(m)*fz)
        call update_stress_ew2p(jatm,3,3,-0.5d0*zdf(m)*fz)
        call update_stress_ew2p(iatm,1,1,-0.5d0*xdf(m)*fx)
        call update_stress_ew2p(iatm,1,2,-0.5d0*xdf(m)*fy)
        call update_stress_ew2p(iatm,1,3,-0.5d0*xdf(m)*fz)
        call update_stress_ew2p(iatm,2,1,-0.5d0*xdf(m)*fy)
        call update_stress_ew2p(iatm,2,2,-0.5d0*ydf(m)*fy)
        call update_stress_ew2p(iatm,2,3,-0.5d0*ydf(m)*fz)
        call update_stress_ew2p(iatm,3,1,-0.5d0*xdf(m)*fz)
        call update_stress_ew2p(iatm,3,2,-0.5d0*ydf(m)*fz)
        call update_stress_ew2p(iatm,3,3,-0.5d0*zdf(m)*fz)
#endif
#endif /* HEAT_CURRENT */

            endif

ccc          endif

        enddo

c
c     load temps back to fxx(iatm) etc

        fxx(iatm) = fi(1)
        fyy(iatm) = fi(2)
        fzz(iatm) = fi(3)

!     for force decompositions

#ifdef TTM_FORCE_DECOMPOSITION
        fx_ttm_per(iatm) = fi_cc(1)
        fy_ttm_per(iatm) = fi_cc(2)
        fz_ttm_per(iatm) = fi_cc(3)

        fx_ttm_ind(iatm) = fi_dc_dd(1)
        fy_ttm_ind(iatm) = fi_dc_dd(2)
        fz_ttm_ind(iatm) = fi_dc_dd(3)
#endif /* TTM_FORCE_DECOMPOSITION */

#ifdef STRESS
c
c     complete stress tensor

        stress(1) = stress(1) + strs1
        stress(2) = stress(2) + strs2
        stress(3) = stress(3) + strs3
        stress(4) = stress(4) + strs2
        stress(5) = stress(5) + strs5
        stress(6) = stress(6) + strs6
        stress(7) = stress(7) + strs3
        stress(8) = stress(8) + strs6
        stress(9) = stress(9) + strs9
#endif
ccc      endif

#ifdef VAMPIR
      call VTEND(100, ierr)
#endif
      return
      end
