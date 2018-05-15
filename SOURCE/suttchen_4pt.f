      subroutine suttchen
     x  (iatm,ik,ntpvdw,engmet,virmet,rmet,dlrpot,ilist,ltype,lstvdw,
     x  ltpvdw,rsqdf,xdf,ydf,zdf,fxx,fyy,fzz,vvv,ggg,stress,rho)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating metal forces atom
c     using sutton chen potential and verlet neighbour list
c     4 point interpolation scheme
c
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - w. smith        june 1995
c     
c     wl
c     2003/05/08 08:45:12
c     1.4
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf),rsqdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      dimension ilist(mxxdf),ltype(mxatms),lstvdw(mxvdw)
      dimension stress(9),rho(mxatms),fi(3)
      dimension ltpvdw(mxvdw)
CDIR$ CACHE_ALIGN fi
#ifdef VAMPIR
      call VTBEGIN(161, ierr)
#endif
c     
c     set cutoff condition for pair forces
      
      rcsq=rmet**2
c     
c     interpolation spacing
      
      rdr = 1.d0/dlrpot
c     
c     initialise potential energy and virial
      
      engmet=0.d0
      virmet=0.d0
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
          ab = ai*(ai-1.d0)*0.5d0 + aj + 0.5d0
        else
          ab = aj*(aj-1.d0)*0.5d0 + ai + 0.5d0
        endif
        
        k0=lstvdw(int(ab))
        k1=ntpvdw+k0+1
        
        if((ltpvdw(k0).ge.100).and.(vvv(1,k0).ne.0.d0)) then 
c     
c     apply truncation of potential
          
          rsq = rsqdf(m)
          
          if(rcsq.gt.rsq)then
            
            rrr = sqrt(rsq)               
            l=int(rrr*rdr)
            ppp=rrr*rdr-dble(l)
            
            
c     calculate interaction energy using 4-point interpolation
            
            vka = vvv(l-1,k0)
            vk0 = vvv(l,k0)
            vk1 = vvv(l+1,k0)
            vk2 = vvv(l+2,k0)
            
            engadd=vk0+
     x        ppp*(-2.d0*vka-3.d0*vk0+6.d0*vk1-vk2+
     x        ppp*(3.d0*(vka-vk0-vk0+vk1)+
     x        ppp*(-vka+vk2+3.d0*(vk0-vk1))))/6.d0
            
            engmet = engadd + engmet
c     
c     calculate forces using 4-point interpolation
            
            gka = ggg(l-1,k0)
            gk0 = ggg(l,k0)
            gk1 = ggg(l+1,k0)
            gk2 = ggg(l+2,k0)
            
            gamma1=gk0+
     x        ppp*(-2.d0*gka-3.d0*gk0+6.d0*gk1-gk2+
     x        ppp*(3.d0*(gka-gk0-gk0+gk1)+
     x        ppp*(-gka+gk2+3.d0*(gk0-gk1))))/6.d0
            
            gka = ggg(l-1,k1)
            gk0 = ggg(l,k1)
            gk1 = ggg(l+1,k1)
            gk2 = ggg(l+2,k1)
            
            gamma2=gk0+
     x        ppp*(-2.d0*gka-3.d0*gk0+6.d0*gk1-gk2+
     x        ppp*(3.d0*(gka-gk0-gk0+gk1)+
     x        ppp*(-gka+gk2+3.d0*(gk0-gk1))))/6.d0
            
            if(ai.gt.aj)then

              gamma = (gamma1-gamma2*(rho(iatm)*vvv(1,k1)+
     x          rho(jatm)*vvv(2,k1)))/rsq

            else

              gamma = (gamma1-gamma2*(rho(iatm)*vvv(2,k1)+
     x          rho(jatm)*vvv(1,k1)))/rsq

            endif
            
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
            
            stress(1) = stress(1)+ xdf(m)*fx
            stress(2) = stress(2)+ xdf(m)*fy
            stress(3) = stress(3)+ xdf(m)*fz
            
            stress(5) = stress(5)+ ydf(m)*fy
            stress(6) = stress(6)+ ydf(m)*fz
            
            stress(9) = stress(9)+ zdf(m)*fz
#endif            
c     
c     calculate virial
            
            virmet=virmet-gamma*rsq
            
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
      
      stress(4)=stress(2)
      stress(7)=stress(3)
      stress(8)=stress(6)
#endif
      
#ifdef VAMPIR
      call VTEND(161, ierr)
#endif
      return
      end



