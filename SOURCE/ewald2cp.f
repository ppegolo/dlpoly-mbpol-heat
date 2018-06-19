#include "assert.h"
! 06 FEB 06 - IUCHI - ESP FOR DMS
! 14 NOV 05 - IUCHI - INTRODUCING FORCE DECOMPOSITION
! 14 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
!
      subroutine ewald2cp
     x  (iatm,ik,engcpe,vircpe,rcut,epsq,ilist,
     x  chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x  dipx,dipy,dipz,alpha,ercp,drewd,polr,polr2,lacs,
     x  efieldkx,efieldky,efieldkz,emux,emuy,emuz,
     x  lthole,athole,nthole,lttm,nttm2,listttm2,
     x  ascc,ascd,lads)
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
!     Last updated: 6 Feb 2006 by S. Iuchi
c     
c***********************************************************************
c     
! from MODULE
#ifdef TTM_FORCE_DECOMPOSITION
      use ttm_forces,  only: fx_ttm_per, fy_ttm_per, fz_ttm_per, 
     $                       fx_ttm_ind, fy_ttm_ind, fz_ttm_ind
#endif /* TTM_FORCE_DECOMPOSITION */
      use ps_type_dms, only: vesp_r, vesp_dc_r, ldms
      
#include "dl_params.inc"

      parameter(g23=1.3541179394264d0)

      logical lthole,lttm,lacs,lads
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ilist(mxxdf)
      dimension chge(mxatms),polr(mxatms),polr2(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension ercp(mxegrd,0:3)
      dimension listttm2(mxatms)
      dimension stress(9)
      dimension fi(3)
#ifdef TTM_FORCE_DECOMPOSITION
      dimension fi_cc(3), fi_dc_dd(3)  ! for force decompositions
#endif /* TTM_FORCE_DECOMPOSITION */

      dimension bbb(0:3),chgchg(0:3),chgdip(0:3),dipdip(0:3)

CDIR$ CACHE_ALIGN fi
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
      
        fi(1) = fxx(iatm)
        fi(2) = fyy(iatm)
        fi(3) = fzz(iatm)
!        
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

! for ESP

          espdum=0.0d0 
          espdcdum=0.0d0
    
c     calculate interatomic distance
            
            rsq=rsqdf(m)
c     
c     apply truncation of potential
            
            if(rcsq.gt.rsq)then

c
c     calculate interatomic distance
              
              rrr = dsqrt(rsq)               

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
c     charge field

              sdc=screendc

              eij_ki = r4pie0*chge(jatm)*(chgdip(1)+sdc*r3i)
              eij_kj = r4pie0*chge(iatm)*(chgdip(1)+sdc*r3i)

c
c     here is the all dipole smeared option, i.e. lads=.true.

              if (lthole .and. lads) then

                bb=(polr(iatm)*polr(jatm))**(1.d0/6.d0)

                assert(bb.gt.1.d-6)

                if (nthole.eq.3) then

                  arob3=ascd*(rrr/bb)**3
                  dum1=dexp(-arob3)

                else if (nthole.eq.4) then

                  arob4=ascd*(rrr/bb)**4
                  dum1=dexp(-arob4)

                else

                  arob=ascd*rrr/bb
                  arob2=arob*arob
                  dum1=dexp(-arob)*(arob2/2.d0+arob+1.d0)

                endif

                eij_ki = eij_ki-r4pie0*chge(jatm)*dum1*r3i
                eij_kj = eij_kj-r4pie0*chge(iatm)*dum1*r3i

              endif

              efieldkx(iatm) = efieldkx(iatm) + eij_ki*xdf(m)
              efieldky(iatm) = efieldky(iatm) + eij_ki*ydf(m)
              efieldkz(iatm) = efieldkz(iatm) + eij_ki*zdf(m)

              efieldkx(jatm) = efieldkx(jatm) - eij_kj*xdf(m)
              efieldky(jatm) = efieldky(jatm) - eij_kj*ydf(m)
              efieldkz(jatm) = efieldkz(jatm) - eij_kj*zdf(m)
c
c     dipole field

              if(polr2(iatm).gt.1.d-6 .and. polr2(jatm).gt.1.d-6) then

                sdd=screendd

                rirj1 = xdf(m)*xdf(m)
                dtens1 = r4pie0*((3.d0*rirj1*r5i-r3i)*sdd
     x                           + rirj1*dipdip(2) - dipdip(1))
                rirj2 = xdf(m)*ydf(m)
                dtens2 = r4pie0*((3.d0*rirj2*r5i)*sdd
     x                           + rirj2*dipdip(2))
                rirj3 = xdf(m)*zdf(m)
                dtens3 = r4pie0*((3.d0*rirj3*r5i)*sdd
     x                           + rirj3*dipdip(2))
                rirj4 = ydf(m)*ydf(m)
                dtens4 = r4pie0*((3.d0*rirj4*r5i-r3i)*sdd
     x                           + rirj4*dipdip(2) - dipdip(1))
                rirj5 = ydf(m)*zdf(m)
                dtens5 = r4pie0*((3.d0*rirj5*r5i)*sdd
     x                           + rirj5*dipdip(2))
                rirj6 = zdf(m)*zdf(m)
                dtens6 = r4pie0*((3.d0*rirj6*r5i-r3i)*sdd
     x                           + rirj6*dipdip(2) - dipdip(1))
c
c      here is the all dipole smeared option, i.e. lads=.true.

                if (lthole .and. lads) then

                  bb=(polr(iatm)*polr(jatm))**(1.d0/6.d0)

                  if (nthole.eq.3) then
                    arob3=athole*(rrr/bb)**3
                    dum1=dexp(-arob3)
                    dum2=(1.d0+arob3)*dum1
                  else if (nthole.eq.4) then
                    arob4=athole*(rrr/bb)**4
                    dum1=dexp(-arob4)
                    dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
                  else
                    arob=athole*rrr/bb
                    arob2=arob*arob
                    arob3=arob2*arob
                    dum1=dexp(-arob)*(arob2/2.d0+arob+1.d0)
                    dum2=dum1+dexp(-arob)*arob3/6.d0
                  endif

                  dtens1=dtens1-(3.d0*rirj1*r5i*dum2-dum1*r3i)*r4pie0
                  dtens2=dtens2-3.d0*rirj2*r5i*dum2*r4pie0
                  dtens3=dtens3-3.d0*rirj3*r5i*dum2*r4pie0
                  dtens4=dtens4-(3.d0*rirj4*r5i*dum2-dum1*r3i)*r4pie0
                  dtens5=dtens5-3.d0*rirj5*r5i*dum2*r4pie0
                  dtens6=dtens6-(3.d0*rirj6*r5i*dum2-dum1*r3i)*r4pie0

                endif

                dmui1=dipx(iatm)
                dmui2=dipy(iatm)
                dmui3=dipz(iatm)
                dmuj1=dipx(jatm)
                dmuj2=dipy(jatm)
                dmuj3=dipz(jatm)

                emux(iatm)=emux(iatm)+dtens1*dmuj1
     x                         +dtens2*dmuj2
     x                         +dtens3*dmuj3
                emuy(iatm)=emuy(iatm)+dtens2*dmuj1
     x                         +dtens4*dmuj2
     x                         +dtens5*dmuj3
                emuz(iatm)=emuz(iatm)+dtens3*dmuj1
     x                         +dtens5*dmuj2
     x                         +dtens6*dmuj3

                emux(jatm)=emux(jatm)+dtens1*dmui1
     x                         +dtens2*dmui2
     x                         +dtens3*dmui3
                emuy(jatm)=emuy(jatm)+dtens2*dmui1
     x                         +dtens4*dmui2
     x                         +dtens5*dmui3
                emuz(jatm)=emuz(jatm)+dtens3*dmui1
     x                         +dtens5*dmui2
     x                         +dtens6*dmui3

              endif
c
c     charge-charge interactions

              cicj=chge(iatm)*chge(jatm)

              edum = cicj*(chgchg(0)+screencc/rrr)*r4pie0
              engcpe=engcpe+edum

              espdum = (chgchg(0)+screencc/rrr)*r4pie0  ! ESP

              ffac = cicj*(screencc*r3i-dscreencc/rsqi+
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

                assert(bb.gt.1.d-6)

                if (nthole.eq.3) then
                  arob3=ascc*(rrr/bb)**3
                  dum1=dexp(-arob3)
                  dum2=dum1*(1.d0-arob3*gammq23(arob3))/rrr*r4pie0
                  engcpe=engcpe-cicj*dum2
                  ffac=ffac-cicj*r3i*dum1*r4pie0
!     for ESP
                  espdum = espdum - dum2
                else if (nthole.eq.4) then
                  arob4=ascc*(rrr/bb)**4
                  dum1=dexp(-arob4)
                  dum2=dum1*(1.d0-arob4*gammq34(arob4))/rrr*r4pie0
                  engcpe=engcpe-cicj*dum2
                  ffac=ffac-cicj*r3i*dum1*r4pie0
!     for ESP
                  espdum = espdum - dum2
                else
                  arob=ascc*rrr/bb
                  arob2=arob*arob
                  dum1=dexp(-arob)*(arob2/2.d0+arob+1.d0)
                  engcpe=engcpe-cicj*dum1/rrr*r4pie0+
     x               0.5d0*cicj*(ascc/bb)*dexp(-arob)*(1+arob)*r4pie0
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

                  arob4=athole*(rrr/bb)**4
                  dum1=dexp(-arob4)
                  dum2=(1.d0+(4.d0/3.d0)*arob4)*dum1
                  dum3=4.d0*athole*(4.d0*arob4-1.d0)*dum1*rsqi/bb**4
                  dum4=4.d0*athole*dum1*rrr/bb**4

                  facmu = facmu-facmu1*dum2-facmu2*dum1+
     x                    rdmui*rdmuj*dum3*rri*r4pie0-
     x                    dpidpj*rsqi*dum4*r4pie0

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

              endif

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
!
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
            endif

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

#ifdef VAMPIR
      call VTEND(100, ierr)
#endif
      return
      end

! original code in ewald1p.f still uses gammq.f

      function gammq23(x)

         implicit none

         real(8) gammq23,x
         real(8) gammcf,gamser

         real(8), parameter :: a = 2.d0/3.d0
         real(8), parameter :: g23 = 1.3541179394264d0

         if(x.lt.1.d0+a)then
            call gser23(gamser,x)
            gammq23=g23*exp(x-a*log(x))-gamser
         else
            call gcf23(gammcf,x)
            gammq23=gammcf
         endif
         return
      end

      subroutine gcf23(gammcf,x)

         implicit none

         real(8) :: gammcf,x,an,b,c,d,del,h

         real(8), parameter :: a = 2.d0/3.d0
         real(8), parameter :: EPS = 3.d-7
         real(8), parameter :: FPMIN = 1.d-30

         integer, parameter :: ITMAX = 100
         integer :: i

         b=x+1.d0-a
         c=1.d0/FPMIN
         d=1.d0/b
         h=d
         do 11 i=1,ITMAX
            an=-i*(i-a)
            b=b+2
            d=an*d+b
            if(abs(d).lt.FPMIN)d=FPMIN
            c=b+an/c
            if(abs(c).lt.FPMIN)c=FPMIN
            d=1.d0/d
            del=d*c
            h=h*del
            if(abs(del-1.d0).lt.EPS)goto 1
 11         continue
            stop
 1       gammcf=h
         return
      end

      subroutine gser23(gamser,x)

         implicit none

         real(8) :: gamser,x,ap,del,sum
         real(8), parameter :: a = 2.d0/3.d0
         real(8), parameter :: EPS = 3.d-7

         integer, parameter :: ITMAX = 800
         integer n

         if(x.le.0.) then
            gamser=0.d0
            return
         endif
         ap=a
         sum=1.d0/a
         del=sum
         do 11 n=1,ITMAX
            ap=ap+1.d0
            del=del*x/ap
            sum=sum+del
            if(abs(del).lt.abs(sum)*EPS)goto 1
 11         continue
         stop ! 'a too large, ITMAX too small in gser'
 1       gamser=sum
         return
      end

      function gammq34(x)

         implicit none

         real(8) gammq34,x
         real(8) gammcf,gamser

         real(8), parameter :: a = 3.d0/4.d0
         real(8), parameter :: g34 = 1.2254167024651963d0

         if(x.lt.1.d0+a)then
            call gser34(gamser,x)
            gammq34=g34*exp(x-a*log(x))-gamser
         else
            call gcf34(gammcf,x)
            gammq34=gammcf
         endif
         return
      end

      subroutine gcf34(gammcf,x)

         implicit none

         real(8) :: gammcf,x,an,b,c,d,del,h

         real(8), parameter :: a = 3.d0/4.d0
         real(8), parameter :: EPS = 3.d-7
         real(8), parameter :: FPMIN = 1.d-30

         integer, parameter :: ITMAX = 100
         integer :: i

         b=x+1.d0-a
         c=1.d0/FPMIN
         d=1.d0/b
         h=d
         do 11 i=1,ITMAX
            an=-i*(i-a)
            b=b+2
            d=an*d+b
            if(abs(d).lt.FPMIN)d=FPMIN
            c=b+an/c
            if(abs(c).lt.FPMIN)c=FPMIN
            d=1.d0/d
            del=d*c
            h=h*del
            if(abs(del-1.d0).lt.EPS)goto 1
 11         continue
            stop
 1       gammcf=h
         return
      end

      subroutine gser34(gamser,x)

         implicit none

         real(8) :: gamser,x,ap,del,sum
         real(8), parameter :: a = 3.d0/4.d0
         real(8), parameter :: EPS = 3.d-7

         integer, parameter :: ITMAX = 800
         integer n

         if(x.le.0.) then
            gamser=0.d0
            return
         endif
         ap=a
         sum=1.d0/a
         del=sum
         do 11 n=1,ITMAX
            ap=ap+1.d0
            del=del*x/ap
            sum=sum+del
            if(abs(del).lt.abs(sum)*EPS)goto 1
 11         continue
         stop ! 'a too large, ITMAX too small in gser'
 1       gamser=sum
         return
      end
