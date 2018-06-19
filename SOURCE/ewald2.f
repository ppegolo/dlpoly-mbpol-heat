      subroutine ewald2
     x  (iatm,ik,engcpe,vircpe,drewd,rcut,epsq,
     x  ilist,chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,erc,fer,stress)
c
c***********************************************************************
c
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c
c     parallel replicated data version (part 2)
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     T3d optimised. t.forester july 1994
c
c     part 2 - real space terms.
c
c     Tabulated potential in r space
c     3pt interpolation
c
c     t. forester March 1993
c     {stress tensor : t.forester june 1994}
c
c     wl
c     2001/06/12 12:47:27
c     1.4
c     Exp
c
c***********************************************************************
c
#ifdef HEAT_CURRENT
      use heatcurrent, only: update_stress_ew2p, update_energy_ew2p,
     x                       update_forces
#endif

#include "dl_params.inc"

      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ilist(mxxdf)
      dimension chge(mxatms),erc(mxegrd),fer(mxegrd)
      dimension stress(9)
      dimension fi(3)

#ifdef HEAT_CURRENT
      real(8), parameter :: third=1.0d0/3.0d0
      real(8) :: force_tmp(3)
#endif

CDIR$ CACHE_ALIGN fi
#ifdef VAMPIR
      call VTBEGIN(100, ierr)
#endif
c
c     set cutoff condition for pair forces

      rcsq=rcut**2
c
c     reciprocal of interpolation interval

      rdrewd = 1.d0/drewd
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
c     initialise potential energy and virial

      engcpe=0.d0
      vircpe=0.d0

c
c     start of primary loop for forces evaluation

      chgea = chge(iatm)/epsq*r4pie0

      if(chgea.ne.0.d0) then

        fi(1) = fxx(iatm)
        fi(2) = fyy(iatm)
        fi(3) = fzz(iatm)

        do m=1,ik
c
c     atomic index and charge product

          jatm=ilist(m)
          chgprd=chgea*chge(jatm)

c
c     Ignore interaction if product of charges is zero

          if(chgprd.ne.0.d0) then
c
c     calculate interatomic distance

            rsq=rsqdf(m)

c
c     apply truncation of potential

            if(rcsq.gt.rsq)then

              rrr = sqrt(rsq)

              ll = int(rrr*rdrewd)
              l1=ll+1
              l2=ll+2
              ppp = rrr*rdrewd - dble(ll)

c
c     calculate interaction energy using 3-point interpolation

              vk = erc(ll)
              vk1 = erc(l1)
              vk2 = erc(l2)

              t1 = vk + (vk1-vk)*ppp
              t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)

              erfcr = (t1 + (t2-t1)*ppp*0.5d0)*chgprd

              engcpe=engcpe+erfcr
#ifdef HEAT_CURRENT
              call update_energy_ew2p(iatm,0.5d0*erfcr)
              call update_energy_ew2p(jatm,0.5d0*erfcr)
#endif /* HEAT_CURRENT */

c
c     calculate forces using 3pt interpolation

              vk = fer(ll)
              vk1 = fer(l1)
              vk2 = fer(l2)

              t1 = vk + (vk1-vk)*ppp
              t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)

              egamma = (t1 + (t2-t1)*ppp*0.5d0)*chgprd
              vircpe = vircpe - egamma*rsq

              fx=egamma*xdf(m)
              fy=egamma*ydf(m)
              fz=egamma*zdf(m)

              fi(1)=fi(1)+fx
              fi(2)=fi(2)+fy
              fi(3)=fi(3)+fz

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

#ifdef HEAT_CURRENT
              force_tmp = (/-fx,-fy,-fz/)
              call update_forces(jatm,iatm,force_tmp)
              call update_forces(iatm,jatm,-force_tmp)
#endif /*HEAT_CURRENT*/

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

          endif

        enddo

c
c     load temps back to fxx(iatm) etc

        fxx(iatm) = fi(1)
        fyy(iatm) = fi(2)
        fzz(iatm) = fi(3)

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
      endif

#ifdef VAMPIR
      call VTEND(100, ierr)
#endif
      return
      end
