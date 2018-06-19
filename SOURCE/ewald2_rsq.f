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
c     
c     part 2 - real space terms. 
c     
c     Tabulated potential in r-squared space
c     2pt interpolation
c     
c     t. forester march 1993
c     t. forester june 1994
c     t.forester stress tensor may 1995
c     
c     wl
c     2001/06/12 12:49:07
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension ilist(mxxdf)
      dimension chge(mxatms),erc(mxegrd),fer(mxegrd)
      dimension stress(9)
      dimension fi(3)
CDIR$ CACHE_ALIGN fi
#ifdef VAMPIR
      call VTBEGIN(101, ierr)
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

      chgea=chge(iatm)/epsq*r4pie0
      if(chgea.ne.0.d0) then

        fi(1) = fxx(iatm)
        fi(2) = fyy(iatm)
        fi(3) = fzz(iatm)
        
c     
c     start of primary loop for forces evaluation
        
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
              
              ll = int(rsq*rdrewd)
              l1=ll+1
              ppp = rsq*rdrewd - dble(ll)
              
c     
c     calculate interaction energy using 2-point interpolation
              
              vk = erc(ll)
              vk1 = erc(l1)
              
              erfcr = (vk + (vk1-vk)*ppp)*chgprd
              
              engcpe=engcpe+erfcr
              
c     
c     calculate forces using 2pt interpolation
              
              vk = fer(ll)
              vk1 = fer(ll+1)
              
              egamma = (vk + (vk1-vk)*ppp)*chgprd
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
      call VTEND(101, ierr)
#endif
      return
      end
