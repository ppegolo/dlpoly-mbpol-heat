      subroutine srfrce
     x  (iatm,ik,engsrp,virsrp,rcut,dlrpot,ilist,ltype,lstvdw,
     x  ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,ggg,stress)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating short range force and
c     potential energy terms using verlet neighbour list
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     
c     version 3
c     author    - t. forester    june  1993
c     stress tensor added t.forester may 1994
c     
c     wl
c     2000/01/18 14:05:56
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension ilist(mxxdf),ltpvdw(mxvdw)
      dimension ltype(mxatms),lstvdw(mxvdw)
      dimension stress(9)
      dimension fi(3)
CDIR$ CACHE_ALIGN fi
#ifdef VAMPIR
      call VTBEGIN(83, ierr)
#endif
c     
c     set cutoff condition for pair forces

      rcsq=rcut**2
c     
c     interpolation spacing
      
      rdr = 1.d0/dlrpot
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
      
      engsrp=0.d0
      virsrp=0.d0
c     
c     store forces for iatm 
      
      ai = dble(ltype(iatm))
      fi(1) = fxx(iatm)
      fi(2) = fyy(iatm)
      fi(3) = fzz(iatm)
      
c     
c     start of primary loop for forces evaluation
      
      do m=1,ik
        
c     
c     atomic and potential function indices
        
        jatm=ilist(m)
        
        aj = dble(ltype(jatm))
        
        if(ai.gt.aj) then
          ab = ai*(ai-1.d0)*0.5d0 + aj+0.5d0
        else
          ab = aj*(aj-1.d0)*0.5d0 + ai+0.5d0
        endif
        
        k=lstvdw(int(ab))
        
        if((ltpvdw(k).lt.100).and.(vvv(1,k).ne.0.d0)) then 
c     
c     apply truncation of potential
          
          rsq = rsqdf(m)
          
          if(rcsq.gt.rsq)then
            
            rrr = sqrt(rsq)               
            l=int(rrr*rdr)
            ppp=rrr*rdr-dble(l)

c            
c     calculate interaction energy using 3-point interpolation
            
            vk = vvv(l,k)
            vk1 = vvv(l+1,k)
            vk2 = vvv(l+2,k)
            
            t1 = vk + (vk1-vk)*ppp
            t2 = vk1 +(vk2 - vk1)*(ppp - 1.0d0)
            
            engsrp = t1 + (t2-t1)*ppp*0.5d0 + engsrp
c     
c     calculate forces using 3-point interpolation
            
            gk = ggg(l,k)
            gk1 = ggg(l+1,k)
            gk2 = ggg(l+2,k)
            
            t1 = gk + (gk1-gk)*ppp
            t2 = gk1 + (gk2-gk1)*(ppp - 1.0d0)
            
            gamma = (t1 +(t2-t1)*ppp*0.5d0)/rsq
c     
c     calculate forces
            
            fx = gamma*xdf(m)
            fy = gamma*ydf(m)
            fz = gamma*zdf(m)
            
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
c     
c     calculate virial
            
            virsrp=virsrp-gamma*rsq
            
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
#ifdef VAMPIR
      call VTEND(83, ierr)
#endif
      return
      end
