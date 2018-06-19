      subroutine coul3neu
     x  (ik,engcpe,vircpe,epsq,rcut,ilist,jlist,
     x  chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,stress)  
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     reaction field  potential
c     Ref: M Neumann, J Chem Phys, 82, 5633, (1985)
c      
c     neutral group implementation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester february 1995
c     stress tensor - t.forester   feb 1995
c     
c     wl
c     2000/01/18 14:05:33
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ilist(mxxdf),jlist(mxxdf)
      dimension chge(mxatms),stress(9)
#ifdef VAMPIR
      call VTBEGIN(93, ierr)
#endif
c     
c     reaction field terms

      b0 = 2.d0*(epsq - 1.d0)/(2.d0*epsq + 1.d0)
      rfld0 = b0/rcut**3
      rfld1 = (1.d0 + b0*0.5d0)/rcut
      rfld2 = rfld0*0.5d0
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
      
      do m=1,ik
c     
c     atomic index and charge product
         
         iatm=ilist(m)
         jatm=jlist(m)
         chgprd=chge(jatm)*chge(iatm)*r4pie0
         
         if(chgprd.ne.0.d0) then
c     
c     calculate interatomic distance
            
            rsq=rsqdf(m)
            rrr = sqrt(rsq)
c               
c     calculate potential energy 
               
            coul0 = chgprd/rrr
            coul = coul0 + chgprd*(rfld2*rsq - rfld1)
            engcpe = engcpe + coul
c     
c     calculate forces and virial
               
            fcoul = coul0/rsq - chgprd*rfld0
            vircpe = vircpe - fcoul*rsq

            fx = fcoul*xdf(m)
            fy = fcoul*ydf(m)
            fz = fcoul*zdf(m)
               
            fxx(iatm)=fxx(iatm)+fx
            fyy(iatm)=fyy(iatm)+fy
            fzz(iatm)=fzz(iatm)+fz
               
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
#ifdef VAMPIR
      call VTEND(93, ierr)
#endif
      return
      end
