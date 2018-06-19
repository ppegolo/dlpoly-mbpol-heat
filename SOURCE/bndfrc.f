! 10 NOV 05 - IUCHI - INTRODUCING FORCE DECOMPOSITION
! 10 NOV 05 - IUCHI - ADD USE MODULE TTM_FORCES
! 28 OCT 05 - IUCHI - ADD EXPONENTIAL TYPE BOND POTENTIAL: KEYBND=7

!SR : added a temporary fix - SHOULD CHANGE IT LATER

!
      subroutine bndfrc
     x  (idnode,imcon,mxnode,ntbond,engbnd,virbnd,keybnd,listbnd,cell,
     x  fxx,fyy,fzz,prmbnd,xxx,yyy,zzz,xdab,ydab,zdab,stress,buffer,
     $  lttm,wat_monomer_ener)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating chemical bond energy and
c     force terms in molecular dynamics.
c
c     replicated data - blocked  data version
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith        july 1992
c
c     'lennard-jones' bonds added for GROMOS 1-4 interactions
c     - t forester     march 1993
c     keybnd > 0 ==> bond (excluded interaction)
c     keybnd < 0 ==> distance restraint (non-excluded interaction)
c     - t. forester    march 1994
c
c     stress tensor added : t.forester may 1994
c     block data version :  t.forester nov 1994
c
c     wl
c     2001/08/31 11:13:42
c     1.12
c     Exp
!
!     Last updated: 10 Nov 2005 by S. Iuchi
c
c***********************************************************************
! from MODULE
#ifdef TTM_FORCE_DECOMPOSITION
      use ttm_forces, only: fx_ttm_intra, fy_ttm_intra, fz_ttm_intra
#endif /* TTM_FORCE_DECOMPOSITION */

#ifdef HEAT_CURRENT
      use heatcurrent, only: update_stress_bnd, update_energy_bnd,
     x                       update_forces
#endif

#include "dl_params.inc"

      logical safe
      logical lttm ! added for force decompositions
      dimension keybnd(mxtbnd),listbnd(mxbond,3)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms),cell(9)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension xdab(msbad),ydab(msbad),zdab(msbad)
      dimension prmbnd(mxtbnd,mxpbnd),buffer(mxbuff)
      dimension stress(9)

      dimension wat_monomer_ener(*)

      integer index_mol

#ifdef HEAT_CURRENT
      real(8), parameter :: half=0.5d0
      real(8) :: force_tmp(3)
#endif

c
c     define chemical bond potential function and derivative
c     using the parameters in array prmbnd

c
c     harmonic bonds

      vb1(r,k)=0.5d0*prmbnd(k,1)*(r-prmbnd(k,2))**2
      gb1(r,k)=prmbnd(k,1)*(r-prmbnd(k,2))

c
c     morse potential bonds

c     original dlpoly: mors
      vb2(r,k)=prmbnd(k,1)*((1.d0-exp(-prmbnd(k,3)*
     x  (r-prmbnd(k,2))))**2-1.d0)
      gb2(r,k)=2.d0*prmbnd(k,1)*prmbnd(k,3)*(1.d0-
     x  exp(-prmbnd(k,3)*(r-prmbnd(k,2))))*
     x  exp(-prmbnd(k,3)*(r-prmbnd(k,2)))

c     for tip4p/2005 flexible: qmrs
      vb8(r,k)=prmbnd(k,1)*((1.d0-exp(-prmbnd(k,3)*
     x  (r-prmbnd(k,2))))**2)
      gb8(r,k)=2.d0*prmbnd(k,1)*prmbnd(k,3)*(1.d0-
     x  exp(-prmbnd(k,3)*(r-prmbnd(k,2))))*
     x  exp(-prmbnd(k,3)*(r-prmbnd(k,2)))

c
c     12-6 potential

      vb3(r,k)=(prmbnd(k,1)*r**(-6) - prmbnd(k,2))*r**(-6)
      gb3(r,k)=(6.d0*prmbnd(k,2) -12.d0*prmbnd(k,1)/r**6)/r**7

c
c     restrained harmonic: ro = r- prmbnd(k,2)

      vb4(ro,k) = 0.5d0*prmbnd(k,1)*(min(abs(ro),prmbnd(k,3)))**2
     x  + prmbnd(k,1)*prmbnd(k,3)*max(abs(ro)-prmbnd(k,3),0.d0)
      gb4(ro,k) = prmbnd(k,1)*(sign(min(abs(ro),prmbnd(k,3)),ro))
c
c     quartic potential

      vb5(r,k)=0.5d0*prmbnd(k,1)*(r-prmbnd(k,2))**2 +
     x  1.d0/3.d0*prmbnd(k,3)*(r-prmbnd(k,2))**3 +
     x  0.25d0*prmbnd(k,4)*(r-prmbnd(k,2))**4
      gb5(r,k)=prmbnd(k,1)*(r-prmbnd(k,2)) +
     x  prmbnd(k,3)*(r-prmbnd(k,2))**2 +
     x  prmbnd(k,4)*(r-prmbnd(k,2))**3
!
!     exponential potential

      vb7(r,k)=prmbnd(k,1)*exp(-prmbnd(k,2)*r)
      gb7(r,k)=-prmbnd(k,1)*prmbnd(k,2)*exp(-prmbnd(k,2)*r)

#ifdef VAMPIR
      call VTBEGIN(21, ierr)
#endif
c
c     flag for undefined potential

      safe = .true.
c
c     check size of work arrays

      if((ntbond-mxnode+1)/mxnode.gt.msbad) call error(idnode,418)
c
c     block indices

      ibnd1 = (idnode*ntbond)/mxnode + 1
      ibnd2 = ((idnode+1)*ntbond)/mxnode

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
c     calculate atom separation vectors

      ii=0
      do i=ibnd1,ibnd2

        ii=ii+1
c
c     indices of atoms involved

        ia=listbnd(ii,2)
        ib=listbnd(ii,3)

c
c     components of bond vector

        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)

      enddo

c
c     periodic boundary condition

      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
c
c     zero bond energy and virial accumulators

      engbnd=0.d0
      virbnd=0.d0

c
c     loop over all specified chemical bond potentials


      ii=0
      do i=ibnd1,ibnd2

        ii=ii+1
c
c     define components of bond vector

        rrab=0.d0
        rab=sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)
        if(rab.gt.1.d-6)rrab=1.d0/rab
c
c     index of potential function parameters

        kk=listbnd(ii,1)
        keyb = abs(keybnd(kk))
c
c     calculate scalar constant terms

        if(keyb.eq.0) then
c
c     null interaction

          omega = 0.d0
          gamma = 0.d0

        elseif(keyb.eq.1)then
c
c     harmonic potential

          omega=vb1(rab,kk)
          gamma=gb1(rab,kk)*rrab

        else if(keyb.eq.2)then
c
c     morse potential

          omega=vb2(rab,kk)
          gamma=gb2(rab,kk)*rrab

        else if(keyb.eq.3) then
c
c     12-6 potential

          omega=vb3(rab,kk)
          gamma=gb3(rab,kk)*rrab

        elseif(keyb.eq.4) then
c
c     restrained harmonic

          rab1 = rab - prmbnd(kk,2)
          omega= vb4(rab1,kk)
          gamma= gb4(rab1,kk)*rrab

        elseif(keyb.eq.5) then
c
c     quartic potential

          omega = vb5(rab,kk)
          gamma = gb5(rab,kk)*rrab

       elseif(keyb.eq.7) then
c
c     exponential potential

          omega = vb7(rab,kk)
          gamma = gb7(rab,kk)*rrab

        else if(keyb.eq.8)then
c
c     tip4p/2005 type morse potential

          omega=vb8(rab,kk)
          gamma=gb8(rab,kk)*rrab

        else
c
c     undefined potential

          omega = 0.d0
          gamma = 0.d0
          safe = .false.

        endif
c
c     calculate bond energy and virial

        engbnd=engbnd+omega
        virbnd=virbnd+gamma*rab*rab

#ifdef HEAT_CURRENT
        call update_energy_bnd(listbnd(ii,2),0.5d0*omega)
        call update_energy_bnd(listbnd(ii,3),0.5d0*omega)
#endif

#ifdef OPT

!SR  *FIX IT LATER*  right now, it assumes water is only the molecule that contains bonds
! water has three bonds - FIELD file
        index_mol=1+(i-1)/3
        wat_monomer_ener(index_mol)=wat_monomer_ener(index_mol)+omega

#endif /* OPT */

c
c     indices of atoms involved

        ia=listbnd(ii,2)
        ib=listbnd(ii,3)
c
c     calculate forces

        fx = -gamma*xdab(ii)
        fy = -gamma*ydab(ii)
        fz = -gamma*zdab(ii)

        fxx(ia) = fxx(ia) + fx
        fyy(ia) = fyy(ia) + fy
        fzz(ia) = fzz(ia) + fz

        fxx(ib) = fxx(ib) - fx
        fyy(ib) = fyy(ib) - fy
        fzz(ib) = fzz(ib) - fz
#ifdef HEAT_CURRENT
        force_tmp = (/ fx,fy,fz /)
        call update_forces(ia,ib,force_tmp)
        call update_forces(ib,ia,-force_tmp)
#endif /*HEAT_CURRENT*/
!
!     for force decompositions
#ifdef TTM_FORCE_DECOMPOSITION
        if( lttm ) then

           fx_ttm_intra(ia) = fx_ttm_intra(ia) + fx
           fy_ttm_intra(ia) = fy_ttm_intra(ia) + fy
           fz_ttm_intra(ia) = fz_ttm_intra(ia) + fz

           fx_ttm_intra(ib) = fx_ttm_intra(ib) - fx
           fy_ttm_intra(ib) = fy_ttm_intra(ib) - fy
           fz_ttm_intra(ib) = fz_ttm_intra(ib) - fz

        endif
#endif /* TTM_FORCE_DECOMPOSITION */

#ifdef STRESS
c
c     calculate stress tensor

        strs1 = strs1 + xdab(ii)*fx
        strs2 = strs2 + xdab(ii)*fy
        strs3 = strs3 + xdab(ii)*fz

        strs5 = strs5 + ydab(ii)*fy
        strs6 = strs6 + ydab(ii)*fz

        strs9 = strs9 + zdab(ii)*fz
#endif

#ifdef HEAT_CURRENT
#ifdef HEAT_STRESS
        call update_stress_bnd(ia,1,1,-half*xdab(ii)*fx)
        call update_stress_bnd(ia,1,2,-half*xdab(ii)*fy)
        call update_stress_bnd(ia,1,3,-half*xdab(ii)*fz)
        call update_stress_bnd(ia,2,1,-half*xdab(ii)*fy)
        call update_stress_bnd(ia,2,2,-half*ydab(ii)*fy)
        call update_stress_bnd(ia,2,3,-half*ydab(ii)*fz)
        call update_stress_bnd(ia,3,1,-half*xdab(ii)*fz)
        call update_stress_bnd(ia,3,2,-half*ydab(ii)*fz)
        call update_stress_bnd(ia,3,3,-half*zdab(ii)*fz)
        call update_stress_bnd(ib,1,1,-half*xdab(ii)*fx)
        call update_stress_bnd(ib,1,2,-half*xdab(ii)*fy)
        call update_stress_bnd(ib,1,3,-half*xdab(ii)*fz)
        call update_stress_bnd(ib,2,1,-half*xdab(ii)*fy)
        call update_stress_bnd(ib,2,2,-half*ydab(ii)*fy)
        call update_stress_bnd(ib,2,3,-half*ydab(ii)*fz)
        call update_stress_bnd(ib,3,1,-half*xdab(ii)*fz)
        call update_stress_bnd(ib,3,2,-half*ydab(ii)*fz)
        call update_stress_bnd(ib,3,3,-half*zdab(ii)*fz)
#endif
#endif
      enddo
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
c
c     check for undefined potentials

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,444)
c
c     sum contributions to potential and virial

      if(mxnode.gt.1) then

        buffer(3)=engbnd
        buffer(4)=virbnd

        call gdsum(buffer(3),2,buffer(1))

        engbnd=buffer(3)
        virbnd=buffer(4)

      endif


#ifdef OPT

      index_mol=ntbond/3

c SR added for wat_monomer_ener

      if(mxnode.gt.1) then
        call gdsum(wat_monomer_ener,index_mol,buffer)
      endif

#endif /* OPT */

#ifdef VAMPIR
      call VTEND(21, ierr)
#endif
      return
      end
