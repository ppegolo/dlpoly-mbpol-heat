      subroutine srfrceneu
     x  (ik,engsrp,virsrp,dlrpot,rcut,ilist,jlist,ltype,
     x  lstvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,ggg,stress)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating short range force and
c     potential energy terms using verlet neighbour list
c     neutral group implementation
c     r-squared interpolation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory  1995
c     author    - t. forester       nov 1995
c     
c     neutral groups
c     author    - t. forester    march  1994
c     
c     wl
c     2000/01/18 14:05:57
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension ilist(mxxdf),jlist(mxxdf)
      dimension ltype(mxatms),lstvdw(mxvdw)
      dimension stress(9)
#ifdef VAMPIR
      call VTBEGIN(116, ierr)
#endif
c     
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
c     
c     reciprocal of interpolation spacing

      rdlpot = 1.d0/dlrpot
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
c     start of primary loop for forces evaluation
      
      do m=1,ik
c     
c     atomic and potential function indices
        
        iatm=ilist(m)
        ai = ltype(iatm)
        jatm=jlist(m)
        aj = ltype(jatm)

        if(ai.gt.aj) then
          ak = (ai*(ai-1.d0)*0.5d0 + aj + 0.5d0)
        else
          ak = (aj*(aj-1.d0)*0.5d0 + ai + 0.5d0)
        endif
        k=lstvdw(int(ak))

        if(vvv(1,k).ne.0.d0) then 

          rsq = rsqdf(m)

          if(rsq.lt.rcsq) then
c     
c     determine interpolation panel for force arrays
            
            l=int(rsq*rdlpot)
            ppp=rsq*rdlpot-dble(l)
c     
c     calculate interaction energy using 2-point interpolation
            
            vk = vvv(l,k)
            vk1 = vvv(l+1,k)
            
            engsrp = vk + (vk1-vk)*ppp + engsrp
c     
c     calculate forces using 2-point interpolation
            
            gk = ggg(l,k)
            gk1 = ggg(l+1,k)
            
            gamma = (gk + (gk1-gk)*ppp)/rsq
            fx = gamma*xdf(m)
            fy = gamma*ydf(m)
            fz = gamma*zdf(m)
            
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
c     
c     calculate virial
            
            virsrp=virsrp-gamma*rsq
            
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
      call VTEND(116, ierr)
#endif
      return
      end
