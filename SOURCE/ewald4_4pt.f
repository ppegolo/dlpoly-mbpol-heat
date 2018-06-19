      subroutine ewald4
     x  (iatm,ik,engcpe,vircpe,engcpl,vircpl,drewd,
     x  rcut,epsq,ilist,chge,rsqdf,xdf,ydf,zdf,
     x  fxx,fyy,fzz,flx,fly,flz,erc,fer,stress,stresl)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c     
c     modified to allow direct calculation of primary (short-range)
c     interactions for multiple-time step corrections
c
c     4pt interpolation scheme
c
c     primary neighbours are taken out of the Ewald sum
c     electrostatics are evaluated directly instead
c     
c     parallel replicated data version (part 2)
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     
c     modified  - t. forester february 1993
c     
c     part 2 - real space terms
c     
c     wl
c     2001/06/12 12:51:28
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension flx(mxatms),fly(mxatms),flz(mxatms)
      dimension ilist(mxxdf)
      dimension chge(mxatms)
      dimension erc(mxegrd),fer(mxegrd)
      dimension stresl(9),stress(9)
      dimension fi(3),fli(3)
CDIR$ CACHE_ALIGN fi
CDIR$ CACHE_ALIGN fli
      
#ifdef VAMPIR
      call VTBEGIN(107, ierr)
#endif
c     
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
c     
c     reciprocal of interpolation interval
      
      rdrewd = 1.d0/drewd
c     
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      engcpl=0.d0
      vircpl=0.d0

#ifdef STRESS
c
c     initialise stress tensor accumulators
      strs1 = 0.d0
      strs2 = 0.d0
      strs3 = 0.d0
      strs5 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0
      strl1 = 0.d0
      strl2 = 0.d0
      strl3 = 0.d0
      strl5 = 0.d0
      strl6 = 0.d0
      strl9 = 0.d0
#endif

c     
c     start of primary loop for forces evaluation
      
      chgea = chge(iatm)/epsq*r4pie0
      if(chgea.ne.0.d0) then
c     
c     temporary arrays for cache aligning
        
        fi(1) = fxx(iatm)
        fi(2) = fyy(iatm)
        fi(3) = fzz(iatm)
        
        fli(1) = flx(iatm)
        fli(2) = fly(iatm)
        fli(3) = flz(iatm)

        do m=1,ik
c     
c     atomic index and charge product
          
          jatm=ilist(m)
          chgprd=chgea*chge(jatm)
          
          if(chgprd.ne.0.d0) then
c     
c     calculate interatomic distance
            
            rsq=rsqdf(m)
            
            if(rcsq.gt.rsq)then
c     
c     coulombic energy
              
              rrr = sqrt(rsq)
              coul = chgprd/rrr
c     
c     calculate Ewald term using 4-point interpolation
              
              ll = int(rrr*rdrewd)
              ppp = rrr*rdrewd - dble(ll)
              
              vka = erc(ll-1)
              vk0 = erc(ll)
              vk1 = erc(ll+1)
              vk2 = erc(ll+2)
              
              engadd=vk0+
     x          ppp*(-2.d0*vka-3.d0*vk0+6.d0*vk1-vk2+
     x          ppp*(3.d0*(vka-vk0-vk0+vk1)+
     x          ppp*(-vka+vk2+3.d0*(vk0-vk1))))/6.d0

              erfcr = engadd*chgprd
c     
c     sum contributions to the totals 
              
              engcpl = engcpl+erfcr-coul
              engcpe = engcpe + coul
              
c     
c     calculate coulombic forces
              
              fcoul = coul/rsq
c     
c     calculate Ewald forces using 4pt interpolation
              
              vka = fer(ll-1)
              vk0 = fer(ll)
              vk1 = fer(ll+1)
              vk2 = fer(ll+2)
              
              engadd=vk0+
     x          ppp*(-2.d0*vka-3.d0*vk0+6.d0*vk1-vk2+
     x          ppp*(3.d0*(vka-vk0-vk0+vk1)+
     x          ppp*(-vka+vk2+3.d0*(vk0-vk1))))/6.d0
              
              egamma = engadd*chgprd
              egamma = egamma -fcoul
              vircpl = vircpl - egamma*rsq
c     
c     add in contributions to the long-range force

              fx = egamma*xdf(m)
              fy = egamma*ydf(m)
              fz = egamma*zdf(m)
              
              fli(1) = fli(1) + fx
              fli(2) = fli(2) + fy
              fli(3) = fli(3) + fz
              
              flx(jatm)=flx(jatm)-fx
              fly(jatm)=fly(jatm)-fy
              flz(jatm)=flz(jatm)-fz
#ifdef STRESS
c     
c     calculate stress tensor
              
              strl1 = strl1 + xdf(m)*fx
              strl2 = strl2 + xdf(m)*fy
              strl3 = strl3 + xdf(m)*fz
              
              strl5 = strl5 + ydf(m)*fy
              strl6 = strl6 + ydf(m)*fz
              
              strl9 = strl9 + zdf(m)*fz
#endif

c     
c     add in contributions to instantaneous force

              vircpe = vircpe - fcoul*rsq
              fx = fcoul*xdf(m)
              fy = fcoul*ydf(m)
              fz = fcoul*zdf(m)
              
              fi(1) = fi(1) + fx
              fi(2) = fi(2) + fy
              fi(3) = fi(3) + fz
              
              fxx(jatm) = fxx(jatm) - fx
              fyy(jatm) = fyy(jatm) - fy
              fzz(jatm) = fzz(jatm) - fz
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

          endif
          
        enddo
c     
c     copy back temporaries

        fxx(iatm) = fi(1)
        fyy(iatm) = fi(2)
        fzz(iatm) = fi(3)

        flx(iatm) = fli(1)
        fly(iatm) = fli(2)
        flz(iatm) = fli(3)

#ifdef STRESS
c     
c     complete stress tensor
        
        stresl(1) = stresl(1) + strl1
        stresl(2) = stresl(2) + strl2
        stresl(3) = stresl(3) + strl3
        stresl(4) = stresl(4) + strl2
        stresl(5) = stresl(5) + strl5
        stresl(6) = stresl(6) + strl6
        stresl(7) = stresl(7) + strl3
        stresl(8) = stresl(8) + strl6
        stresl(9) = stresl(9) + strl9
        
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
      call VTEND(107, ierr)
#endif
      return
      end
