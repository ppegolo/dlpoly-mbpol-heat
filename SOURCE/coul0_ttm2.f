      subroutine coul0_ttm2
     x  (iatm,ik,engcpe,vircpe,rcut,epsq,ilist,
     x  chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)  
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     1/r potential, no truncation or damping
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester february 1993
c     stress tensor - t.forester may 1994
c     
c     wl
c     2000/01/18 14:05:32
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ilist(mxxdf)
      dimension chge(mxatms),stress(9)
      dimension fi(3)

CDIR$ CACHE_ALIGN fi
      
#ifdef VAMPIR
      call VTBEGIN(84, ierr)
#endif
c     
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
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
          
          if(chgprd.ne.0.d0) then
c     
c     calculate interatomic distance
            
            rsq=rsqdf(m)
c     
c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr = sqrt(rsq)
c              
c     calculate potential energy and virial
              
              coul = chgprd/rrr
              engcpe = engcpe + coul
              
c     
c     calculate forces
              
              fcoul = coul/rsq
              fx = fcoul*xdf(m)
              fy = fcoul*ydf(m)
              fz = fcoul*zdf(m)
              
              fi(1)=fi(1)+fx
              fi(2)=fi(2)+fy
              fi(3)=fi(3)+fz
              
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
c     load temps back to fxx(iatm) etc
        
        fxx(iatm) = fi(1)
        fyy(iatm) = fi(2)
        fzz(iatm) = fi(3)
c     
c     virial
        
        vircpe = -engcpe
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
      call VTEND(84, ierr)
#endif
      return
      end
