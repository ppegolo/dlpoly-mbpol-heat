      subroutine ewald2neu
     x  (ik,engcpe,vircpe,drewd,rcut,epsq,ilist,jlist,
     x  chge,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,erc,fer,stress)
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
c     charge groups version  t.forester june 1996
c     
c     part 2 - real space terms. 
c     with charge group implementation
c     
c     Tabulated potential in r space
c     3pt interpolation
c     
c     t. forester March 1993
c     {stress tensor : t.forester june 1994}
c     
c     wl
c     2001/06/12 12:47:58
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ilist(mxxdf),jlist(mxxdf)
      dimension chge(mxatms),erc(mxegrd),fer(mxegrd)
      dimension stress(9)
#ifdef VAMPIR
      call VTBEGIN(137, ierr)
#endif
c     
c     set cutoff condition for pair forces
      
      rcsq=(drewd*dble(mxegrd-4))**2
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
      
      reps = r4pie0/epsq
      do m =1,ik
c     
c     atomic index and charge product
        
        iatm=ilist(m)
        jatm=jlist(m)
        chgprd=chge(jatm)*chge(iatm)*reps
        
        if(chgprd.ne.0.d0) then
          
          rsq = rsqdf(m)
c     
c     apply truncation of potential
          
          if(rsq.lt.rcsq) then

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
            
            fxx(iatm)=fxx(iatm)+fx
            fyy(iatm)=fyy(iatm)+fy
            fzz(iatm)=fzz(iatm)+fz
            
            fxx(jatm)=fxx(jatm)-fx
            fyy(jatm)=fyy(jatm)-fy
            fzz(jatm)=fzz(jatm)-fz
            
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
      call VTEND(137, ierr)
#endif
      return
      end

