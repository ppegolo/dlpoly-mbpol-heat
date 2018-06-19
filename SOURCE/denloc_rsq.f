      subroutine denloc
     x  (iatm,ik,ntpvdw,rvdw,dlrpot,ilist,ltype,lstvdw,ltpvdw,rsqdf,
     x  rho,vvv)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating local atomic density
c     for sutton-chen metal potential 
c     r-squared interpolation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - w. smith       june  1995
c     
c     
c     wl
c     2003/05/08 08:45:10
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension vvv(mxgrid,mxvdw),rsqdf(mxxdf),rho(mxatms)
      dimension ilist(mxxdf),ltype(mxatms),lstvdw(mxvdw),ltpvdw(mxvdw)

#ifdef VAMPIR
      call VTBEGIN(111, ierr)
#endif
c     
c     set cutoff condition for density calculation
      
      rcsq=rvdw**2

c     
c     interpolation spacing

      rdr = 1.d0/dlrpot

c     
c     start of primary loop for density
      
      ai=dble(ltype(iatm))

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
        
        if((ltpvdw(k0).ge.100).and.(vvv(1,k1).ne.0.d0)) then 
c     
c     apply truncation of potential
          
          rsq = rsqdf(m)
          
          if(rcsq.gt.rsq)then
            
            l=int(rsq*rdr)
            ppp=rsq*rdr-dble(l)
            
            
c     calculate density using 2-point interpolation
            
            vk0 = vvv(l,k1)
            vk1 = vvv(l+1,k1)
            
            density = vk0 + (vk1-vk0)*ppp

            if(ai.gt.aj)then

              rho(iatm) = rho(iatm) + density*vvv(1,k1)
              rho(jatm) = rho(jatm) + density*vvv(2,k1)
              
            else

              rho(iatm) = rho(iatm) + density*vvv(2,k1)
              rho(jatm) = rho(jatm) + density*vvv(1,k1)

            endif

          endif
          
        endif
        
      enddo

#ifdef VAMPIR
      call VTEND(111, ierr)
#endif
      return
      end



