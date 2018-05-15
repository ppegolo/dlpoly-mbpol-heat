! 10 NOV 05 - IUCHI - INTRODUCING FORCE DECOMPOSITION
! 10 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
! 31 OCT 05 - IUCHI - DIMENSION FDUM(3,2)
! 31 OCT 05 - IUCHI - CALL PS_TYPE_COUPLE
! 31 OCT 05 - IUCHI - ADD P-S TYPE BOND-ANGLE COUPLE TERM: KEYANG=9
!
!SR : added a temporary fix - mbpol ceiling - SHOULD CHANGE IT LATER

      subroutine angfrc
     x  (idnode,imcon,mxnode,ntangl,engang,virang,keyang,listang,
     x  cell,fxx,fyy,fzz,prmang,xxx,yyy,zzz,xdab,ydab,zdab,xdbc,
     x  ydbc,zdbc,stress,buffer,lttm,wat_monomer_ener)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating bond angle energy and
c     force terms in molecular dynamics.
c
c     replicated data - blocked version
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith         may 1992
c     modified  - t. forester      feb 1993
c     modified  - t.forester       nov 1994 : block data
c     stress tensor -  t.forester  may 1995
c
c     wl
c     2001/05/30 12:39:59
c     1.6
c     Exp
!
!     Last updated: 10 Nov 2005 by S. Iuchi
c
c***********************************************************************
c
! from MODULE
#ifdef TTM_FORCE_DECOMPOSITION
      use ttm_forces, only: fx_ttm_intra, fy_ttm_intra, fz_ttm_intra
#endif /* TTM_FORCE_DECOMPOSITION */

#ifdef HEAT_CURRENT
      use heatcurrent, only: update_stress_ang, update_energy_ang
#endif

#include "dl_params.inc"

      logical safe
      logical lttm ! added for force decompositions
      dimension keyang(mxtang),listang(mxangl,4),prmang(mxtang,mxpang)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdab(msbad),ydab(msbad),zdab(msbad)
      dimension xdbc(msbad),ydbc(msbad),zdbc(msbad)
      dimension cell(9),stress(9),buffer(mxbuff)
      dimension fdum(3,2)
      dimension wat_monomer_ener(*)

      real(kind=8) :: tmp_wat_monomer_ener(ntangl)
      integer                  :: nowater

#ifdef HEAT_CURRENT
      real(8), parameter :: third=1.0d0/3.0d0
      real(8), parameter :: half=1.0d0/2.0d0
#endif


c
c     define angular potential function and derivative
c     using the parameters in array prmang
c
c     Harmonic potential

      vt1(theta,k)=0.5d0*prmang(k,1)*(theta-prmang(k,2))**2
      gt1(theta,k)=prmang(k,1)*(theta-prmang(k,2))
c
c     Quartic potential

      vt2(theta,k)=0.5d0*prmang(k,1)*(theta-prmang(k,2))**2 +
     x  1.d0/3.d0*prmang(k,3)*(theta-prmang(k,2))**3 +
     x  0.25d0*prmang(k,4)*(theta-prmang(k,2))**4
      gt2(theta,k)=prmang(k,1)*(theta-prmang(k,2)) +
     x  prmang(k,3)*(theta-prmang(k,2))**2 +
     x  prmang(k,4)*(theta-prmang(k,2))**3

c
c     truncation and screening functions

      sw1(x,y,a)=exp(-(x**8+y**8)/a**8)
      sw2(x,y,a,b)=exp(-(x/a+y/b))
c
c     truncated harmonic valence angle potential

      vt3(t,x,y,a,b,c)=0.5d0*a*(t-b)**2*sw1(x,y,c)
      gt3(t,x,y,a,b,c)=a*(t-b)*sw1(x,y,c)
c
c     screened harmonic valence angle potential

      vt4(t,x,y,a,b,c,d)=0.5d0*a*(t-b)**2*sw2(x,y,c,d)
      gt4(t,x,y,a,b,c,d)=a*(t-b)*sw2(x,y,c,d)
c
c     screened vessal potential type 1

      vt5(t,x,y,a,b,c,d) = a/(8.d0*(b-pi)**2)*
     x  (((b-pi)**2 - (t-pi)**2)**2)*sw2(x,y,c,d)
      gt5(t,x,y,a,b,c,d) = a/(2.d0*(b-pi)**2)*
     x  ((b-pi)**2 - (t-pi)**2)*(t-pi)*sw2(x,y,c,d)
c
c     truncated vessal potential type 2

      vt6(t,x,y,a,b,c,d)=a*(t**c*(t-b)**2*(t+b-2.d0*pi)**2-
     x  0.5d0*c*pi**(c-1.d0)*(t-b)**2*(pi-b)**3)*sw1(x,y,d)
      gt6(t,x,y,a,b,c,d)=a*(t**(c-1.d0)*(t-b)*(t+b-2.d0*pi)*
     x  ((c+4.d0)*t**2-2.d0*pi*(c+2.d0)*t+c*b*(2.d0*pi-b))-
     x  c*pi**(c-1.d0)*(t-b)*(pi-b)**3)*sw1(x,y,d)
c
c     harmonic cosine potential (note cancellation of sint in gt7)

      vt7(theta,k)=0.5d0*prmang(k,1)*(cos(theta)-cos(prmang(k,2)))**2
      gt7(theta,k)=-prmang(k,1)*(cos(theta)-cos(prmang(k,2)))
c
c     ordinary cosine potential

      vt8(theta,k)=prmang(k,1)*(1+cos(prmang(k,3)*theta-prmang(k,2)))
      gt8(theta,k)=-prmang(k,1)*prmang(k,3)*sin(prmang(k,3)*theta-
     x  prmang(k,2))

#ifdef VAMPIR
      call VTBEGIN(22, ierr)
#endif
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
c SR : initialize water energy array

      tmp_wat_monomer_ener=0.d0

c
c     flag for undefined potentials

      safe = .true.
c
c     check size of work arrays

      if((ntangl-mxnode+1)/mxnode.gt.msbad) call error(idnode,419)
c
c     block indices

      iang1 = (idnode*ntangl)/mxnode + 1
      iang2 = ((idnode+1)*ntangl)/mxnode
c
c     calculate atom separation vectors

      ii=0
      do i=iang1,iang2

        ii=ii+1
c
c     indices of atoms involved

        ia=listang(ii,2)
        ib=listang(ii,3)
        ic=listang(ii,4)

c
c     components of first bond vector

        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)

c
c     components of second bond vector

        xdbc(ii)=xxx(ic)-xxx(ib)
        ydbc(ii)=yyy(ic)-yyy(ib)
        zdbc(ii)=zzz(ic)-zzz(ib)

      enddo

c
c     periodic boundary condition

      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      call images(imcon,0,1,ii,cell,xdbc,ydbc,zdbc)
c
c     zero angle energy accumulator

      engang=0.d0
      virang=0.d0

c
c     loop over all specified angle potentials

      ii=0
      do i=iang1,iang2

        ii=ii+1
c
c     define components of first bond vector

        rab = sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)
        rrab = 1.d0/rab

        xab=xdab(ii)*rrab
        yab=ydab(ii)*rrab
        zab=zdab(ii)*rrab

c
c     define components of second bond vector

        rbc = sqrt(xdbc(ii)**2+ydbc(ii)**2+zdbc(ii)**2)
        rrbc = 1.d0/rbc

        xbc=xdbc(ii)*rrbc
        ybc=ydbc(ii)*rrbc
        zbc=zdbc(ii)*rrbc

c
c     index of potential function parameters

        kk=listang(ii,1)

c
c     determine bond angle and calculate potential energy

        cost=(xab*xbc+yab*ybc+zab*zbc)
        theta=acos(cost)
        sint = max(1.d-8,sqrt(1.d0-cost**2))

        keya = abs(keyang(kk))

        if(keya.eq.1) then
c
c     Harmonic potential

          pterm = vt1(theta,kk)
          vterm = 0.d0
          gamma=  gt1(theta,kk)/sint
          gamsa= 0.d0
          gamsc= 0.d0

        elseif (keya.eq.2) then
c
c     Quartic potential

          pterm = vt2(theta,kk)
          vterm = 0.d0
          gamma = gt2(theta,kk)/sint
          gamsa= 0.d0
          gamsc= 0.d0

        elseif(keya.eq.3)then
c
c     truncated Harmonic potential

          pterm=vt3(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3))
          vterm=-8.d0*pterm*(rab**8+rbc**8)/
     x      prmang(kk,3)**8
          gamma=gt3(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3))/sint
          gamsa=(8.d0*pterm/prmang(kk,3)**8)*rab**7
          gamsc=(8.d0*pterm/prmang(kk,3)**8)*rbc**7

        elseif(keya.eq.4)then
c
c     screened Harmonic potential

          pterm=vt4(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3),
     x      prmang(kk,4))
          vterm=-pterm*(rab/prmang(kk,3)+
     x      rbc/prmang(kk,4))
          gamma=gt4(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3),
     x      prmang(kk,4))/sint
          gamsa=(pterm/prmang(kk,3))
          gamsc=(pterm/prmang(kk,4))

        elseif(keya.eq.5)then
c
c     screened vessal potential (type 1)

          pterm=vt5(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3),
     x      prmang(kk,4))
          vterm=-pterm*(rab/prmang(kk,3)+
     x      rbc/prmang(kk,4))
          gamma=gt5(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3),
     x      prmang(kk,4))/sint
          gamsa=(pterm/prmang(kk,3))
          gamsc=(pterm/prmang(kk,4))

        elseif(keya.eq.6)then
c
c     Truncated Vessal potential (type 2)

          pterm=vt6(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3),
     x      prmang(kk,4))
          vterm=-8.d0*pterm*(rab**8+rbc**8)/
     x      prmang(kk,4)**8
          gamma=gt6(theta,rab,rbc,prmang(kk,1),
     x      prmang(kk,2),prmang(kk,3),
     x      prmang(kk,4))/sint
          gamsa=(8.d0*pterm/prmang(kk,4)**8)*rab**7
          gamsc=(8.d0*pterm/prmang(kk,4)**8)*rbc**7

        elseif(keya.eq.7)then
c
c     harmonic cosine potential

          pterm = vt7(theta,kk)
          vterm = 0.d0
          gamma = gt7(theta,kk)
          gamsa = 0.d0
          gamsc = 0.d0

        elseif(keya.eq.8)then
c
c     ordinary cosine potential

          pterm = vt8(theta,kk)
          vterm = 0.d0
          gamma = gt8(theta,kk)/sint
          gamsa = 0.d0
          gamsc = 0.d0

        elseif(keya.eq.9)then
c
c     couple terms of Partridge and Schwenke type potential

           call ps_type_couple( rab, xab, yab, zab, rbc, xbc, ybc, zbc,
     $          theta, pterm, fdum, vterm )

           fxa = fdum(1,1)
           fya = fdum(2,1)
           fza = fdum(3,1)

           fxc = fdum(1,2)
           fyc = fdum(2,2)
           fzc = fdum(3,2)

!           write(6,*) 'CHECK', pterm * 0.01d0 * 2.39006d-1 ! kcal/mol

        else
c
c     undefined potential

          safe = .false.
          pterm=0.d0
          vterm=0.d0
          gamma=0.d0
          gamsa = 0.d0
          gamsc = 0.d0

        endif

        engang=engang+pterm

#ifdef HEAT_CURRENT
        ia=listang(ii,2)
        ib=listang(ii,3)
        ic=listang(ii,4)
        call update_energy_ang(ia,pterm/3)
        call update_energy_ang(ib,pterm/3)
        call update_energy_ang(ic,pterm/3)
#endif

#ifdef OPT
!SR  for water monomer mbpol ! SHOULD BE CHANGED LATER

        tmp_wat_monomer_ener(i)=tmp_wat_monomer_ener(i)+pterm

#endif /* OPT */
c
c     indices of atoms involved

        ia=listang(ii,2)
        ib=listang(ii,3)
        ic=listang(ii,4)

        if( keya == 9 ) goto 999  ! for PS type
c
c     calculate atomic forces

        fxa = gamma*(xbc-xab*cost)*rrab+gamsa*xab
        fya = gamma*(ybc-yab*cost)*rrab+gamsa*yab
        fza = gamma*(zbc-zab*cost)*rrab+gamsa*zab

        fxc = gamma*(xab-xbc*cost)*rrbc+gamsc*xbc
        fyc = gamma*(yab-ybc*cost)*rrbc+gamsc*ybc
        fzc = gamma*(zab-zbc*cost)*rrbc+gamsc*zbc

 999    continue                ! for PS type PES

        fxx(ia)=fxx(ia)+fxa
        fyy(ia)=fyy(ia)+fya
        fzz(ia)=fzz(ia)+fza

        fxx(ib)=fxx(ib)-fxa-fxc
        fyy(ib)=fyy(ib)-fya-fyc
        fzz(ib)=fzz(ib)-fza-fzc

        fxx(ic)=fxx(ic)+fxc
        fyy(ic)=fyy(ic)+fyc
        fzz(ic)=fzz(ic)+fzc
!
!     for force decompositions

#ifdef TTM_FORCE_DECOMPOSITION
        if( lttm ) then

           fx_ttm_intra(ia) = fx_ttm_intra(ia) + fxa
           fy_ttm_intra(ia) = fy_ttm_intra(ia) + fya
           fz_ttm_intra(ia) = fz_ttm_intra(ia) + fza

           fx_ttm_intra(ib) = fx_ttm_intra(ib) - fxa - fxc
           fy_ttm_intra(ib) = fy_ttm_intra(ib) - fya - fyc
           fz_ttm_intra(ib) = fz_ttm_intra(ib) - fza - fzc

           fx_ttm_intra(ic) = fx_ttm_intra(ic) + fxc
           fy_ttm_intra(ic) = fy_ttm_intra(ic) + fyc
           fz_ttm_intra(ic) = fz_ttm_intra(ic) + fzc

        end if
#endif /* TTM_FORCE_DECOMPOSITION */
c
c     virial

        virang=virang + vterm

#ifdef STRESS
c
c     calculate stress tensor

        strs1 = strs1 + rab*xab*fxa + rbc*xbc*fxc
        strs2 = strs2 + rab*xab*fya + rbc*xbc*fyc
        strs3 = strs3 + rab*xab*fza + rbc*xbc*fzc
        strs5 = strs5 + rab*yab*fya + rbc*ybc*fyc
        strs6 = strs6 + rab*yab*fza + rbc*ybc*fzc
        strs9 = strs9 + rab*zab*fza + rbc*zbc*fzc
#endif

#ifdef HEAT_CURRENT
        call update_stress_ang(ia,1,1,half*(rab*xab*fxa))
        call update_stress_ang(ia,1,2,half*(rab*xab*fya))
        call update_stress_ang(ia,1,3,half*(rab*xab*fza))
        call update_stress_ang(ia,2,1,half*(rab*yab*fxa))
        call update_stress_ang(ia,2,2,half*(rab*yab*fya))
        call update_stress_ang(ia,2,3,half*(rab*yab*fza))
        call update_stress_ang(ia,3,1,half*(rab*zab*fxa))
        call update_stress_ang(ia,3,2,half*(rab*zab*fya))
        call update_stress_ang(ia,3,3,half*(rab*zab*fza))
        call update_stress_ang(ib,1,1,half*(rab*xab*fxa+rbc*xbc*fxc))
        call update_stress_ang(ib,1,2,half*(rab*xab*fya+rbc*xbc*fyc))
        call update_stress_ang(ib,1,3,half*(rab*xab*fza+rbc*xbc*fzc))
        call update_stress_ang(ib,2,1,half*(rab*yab*fxa+rbc*ybc*fxc))
        call update_stress_ang(ib,2,2,half*(rab*yab*fya+rbc*ybc*fyc))
        call update_stress_ang(ib,2,3,half*(rab*yab*fza+rbc*ybc*fzc))
        call update_stress_ang(ib,3,1,half*(rab*zab*fxa+rbc*zbc*fxc))
        call update_stress_ang(ib,3,2,half*(rab*zab*fya+rbc*zbc*fyc))
        call update_stress_ang(ib,3,3,half*(rab*zab*fza+rbc*zbc*fzc))
        call update_stress_ang(ic,1,1,half*(rbc*xbc*fxc))
        call update_stress_ang(ic,1,2,half*(rbc*xbc*fyc))
        call update_stress_ang(ic,1,3,half*(rbc*xbc*fzc))
        call update_stress_ang(ic,2,1,half*(rbc*ybc*fxc))
        call update_stress_ang(ic,2,2,half*(rbc*ybc*fyc))
        call update_stress_ang(ic,2,3,half*(rbc*ybc*fzc))
        call update_stress_ang(ic,3,1,half*(rbc*zbc*fxc))
        call update_stress_ang(ic,3,2,half*(rbc*zbc*fyc))
        call update_stress_ang(ic,3,3,half*(rbc*zbc*fzc))
#endif

      enddo

c
c     check for undefined potentials

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,440)
c
c     sum up contributions to potential and virial

      if (mxnode.gt.1) then

        buffer(1)=engang
        buffer(2)=virang
        call gdsum(buffer(1),2,buffer(3))
        engang=buffer(1)
        virang=buffer(2)

      endif

#ifdef OPT

c  SR added for wat_monomer_ener

      if (mxnode.gt.1) then
        call gdsum(tmp_wat_monomer_ener,ntangl,buffer)
      endif

       wat_monomer_ener(1:ntangl)=wat_monomer_ener(1:ntangl)+tmp_wat_monomer_ener

#endif  /* OPT */



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
      call VTEND(22, ierr)
#endif

      return
      end
