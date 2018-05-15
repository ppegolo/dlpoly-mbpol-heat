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
c     author    - w. smith march 1992.
c     stress tensor t.forester june 1994
c     
c     wl
c     2000/01/18 14:05:56
c     1.3
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
#ifdef VAMPIR
      call VTBEGIN(114, ierr)
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
      
      engsrp=0.d0
      virsrp=0.d0
c     
c     start of primary loop for forces evaluation
      
      do m=1,ik
c     
c     atomic and potential function indices
        
        jatm=ilist(m)
        ka=max(ltype(iatm),ltype(jatm))
        kb=min(ltype(iatm),ltype(jatm))
        k=lstvdw((ka*(ka-1))/2+kb)
        
        if((ltpvdw(k).lt.100).and.(vvv(1,k).ne.0.d0)) then
c     
c     calculate square of interatomic distance
          
          rsq=rsqdf(m)
c     
c     apply truncation of potential
          
          if(rcsq.gt.rsq)then
c     
c     determine interpolation panel for force arrays
            
            rrr=sqrt(rsq)
            rsqi=1.d0/rsq            
            l=int(rrr/dlrpot)
            ppp=rrr/dlrpot-dble(l)
            la = l-1
            l1 = l+1
            l2 = l+2
            vka = vvv(la,k)
            vk0 = vvv(l,k)
            vk1 = vvv(l1,k)
            vk2 = vvv(l2,k)
c     
c     calculate interaction energy using 4-point interpolation
            
            engadd=vk0+
     x        ppp*(-2.d0*vka-3.d0*vk0+6.d0*vk1-vk2+
     x        ppp*(3.d0*(vka-vk0-vk0+vk1)+
     x        ppp*(-vka+vk2+3.d0*(vk0-vk1))))/6.d0
            engsrp = engsrp + engadd
            
c     calculate forces using 4-point interpolation

            gka = ggg(la,k)
            gk0 = ggg(l,k)
            gk1 = ggg(l1,k)
            gk2 = ggg(l2,k)
            
            gamma=rsqi*(gk0+
     x        ppp*(-2.d0*gka-3.d0*gk0+6.d0*gk1-gk2+
     x        ppp*(3.d0*(gka-gk0-gk0+gk1)+
     x        ppp*(-gka+gk2+3.d0*(gk0-gk1))))/6.d0)
            
            fx = gamma*xdf(m)
            fy = gamma*ydf(m)
            fz = gamma*zdf(m)
            
            fxx(iatm) = fxx(iatm) + fx
            fyy(iatm) = fyy(iatm) + fy
            fzz(iatm) = fzz(iatm) + fz
            
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
      call VTEND(114, ierr)
#endif
      return
      end
