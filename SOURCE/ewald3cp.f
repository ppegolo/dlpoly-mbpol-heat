! 06 FEB 06 - IUCHI - ESP FOR DMS
! 14 NOV 05 - IUCHI - INTRODUCING FORCE DECOMPOSITION
! 14 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
!
      subroutine ewald3cp
     x   (iatm,ilst,engcpe,vircpe,alpha,epsq,nexatm,
     x   lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress,
     x   dipx,dipy,dipz,polr,polr2,athole12,athole13,
     x   nthole,lthole,lttm,nttm2,listttm2,
     x   efieldkx,efieldky,efieldkz,emux,emuy,emuz)
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
!     Last updated: 6 Feb 2006 by S. Iuchi
c     
c***********************************************************************
c     
! from MODULE
#ifdef TTM_FORCE_DECOMPOSITION
      use ttm_forces, only: fx_ttm_per, fy_ttm_per, fz_ttm_per, 
     $                      fx_ttm_ind, fy_ttm_ind, fz_ttm_ind
#endif /* TTM_FORCE_DECOMPOSITION */
      use ps_type_dms, only: vesp_c, vesp_dc_c, ldms
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension chge(mxatms),polr(mxatms),polr2(mxatms)
      dimension efieldkx(mxatms),efieldky(mxatms),efieldkz(mxatms)
      dimension emux(mxatms),emuy(mxatms),emuz(mxatms)
      dimension dipx(mxatms),dipy(mxatms),dipz(mxatms)
      dimension listttm2(mxatms)
      dimension stress(9)

      dimension bbb(0:3),chgchg(0:3),chgdip(0:3),dipdip(0:3)

      logical lthole,lttm

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
c     charge field

        eij_ki = r4pie0*chge(jatm)*chgdip(1)/epsq
        eij_kj = r4pie0*chge(iatm)*chgdip(1)/epsq

        efieldkx(iatm) = efieldkx(iatm) + eij_ki*xdf(m)
        efieldky(iatm) = efieldky(iatm) + eij_ki*ydf(m)
        efieldkz(iatm) = efieldkz(iatm) + eij_ki*zdf(m)

        efieldkx(jatm) = efieldkx(jatm) - eij_kj*xdf(m)
        efieldky(jatm) = efieldky(jatm) - eij_kj*ydf(m)
        efieldkz(jatm) = efieldkz(jatm) - eij_kj*zdf(m)
c
c     add intramolecular dipole-dipole stuff
c
        if(polr2(iatm).gt.1.d-6 .and. polr2(jatm).gt.1.d-6) then

            rirj1 = xdf(m)*xdf(m)
            dtens1 = rirj1*dipdip(2) - dipdip(1)
 
            rirj2 = xdf(m)*ydf(m)
            dtens2 = rirj2*dipdip(2)

            rirj3 = xdf(m)*zdf(m)
            dtens3 = rirj3*dipdip(2)

            rirj4 = ydf(m)*ydf(m)
            dtens4 = rirj4*dipdip(2) - dipdip(1)

            rirj5 = ydf(m)*zdf(m)
            dtens5 = rirj5*dipdip(2)

            rirj6 = zdf(m)*zdf(m)
            dtens6 = rirj6*dipdip(2) - dipdip(1)
c
c    Applequist's point dipole model

            dtens1=dtens1+3.d0*rirj1*r5i-r3i
            dtens2=dtens2+3.d0*rirj2*r5i
            dtens3=dtens3+3.d0*rirj3*r5i
            dtens4=dtens4+3.d0*rirj4*r5i-r3i
            dtens5=dtens5+3.d0*rirj5*r5i
            dtens6=dtens6+3.d0*rirj6*r5i-r3i

            if (lthole) then

                bb=(polr(iatm)*polr(jatm))**(1.d0/6.d0)
               athole=athole13
               if(mod(iatm-1,3).eq.0.or.mod(jatm-1,3).eq.0)
     x             athole=athole12

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

                dtens1=dtens1-3.d0*rirj1*r5i*dum2+dum1*r3i
                dtens2=dtens2-3.d0*rirj2*r5i*dum2
                dtens3=dtens3-3.d0*rirj3*r5i*dum2
                dtens4=dtens4-3.d0*rirj4*r5i*dum2+dum1*r3i
                dtens5=dtens5-3.d0*rirj5*r5i*dum2
                dtens6=dtens6-3.d0*rirj6*r5i*dum2+dum1*r3i

             endif

             dmui1=dipx(iatm)*r4pie0/epsq
             dmui2=dipy(iatm)*r4pie0/epsq
             dmui3=dipz(iatm)*r4pie0/epsq
             dmuj1=dipx(jatm)*r4pie0/epsq
             dmuj2=dipy(jatm)*r4pie0/epsq
             dmuj3=dipz(jatm)*r4pie0/epsq

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


        endif ! dipole-dipole
c
c     remove charge-charge interactions

           cicj=chge(iatm)*chge(jatm)

           engcpe = engcpe+cicj*chgchg(0)*r4pie0/epsq
!
!     for ESP
           espdum   = chgchg(0) * r4pie0 / epsq  
           espdcdum = chgdip(1) * r4pie0 / epsq

           ffac = cicj*chgchg(1)*r4pie0/epsq

           ffac_cc = cicj*chgchg(1)*r4pie0/epsq ! for force decompositions
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
c     add dipole-dipole interactions

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

               athole=athole13
               if(mod(iatm-1,3).eq.0.or.mod(jatm-1,3).eq.0)
     x             athole=athole12

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
     x                    rdmui*rdmuj*dum3*rri*r4pie0/epsq-
     x                    dpidpj*rsqi*dum4*r4pie0/epsq
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

           endif ! dipole-dipole

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

