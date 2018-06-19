      subroutine ewald3_neu
     x   (ilst,engcpe,vircpe,alpha,epsq,nexatm,
     x   lexatm,chge,xdf,ydf,zdf,fxx,fyy,fzz,stress)
c     
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c     based on charge group decomposition strategy
c     
c     parallel replicated data version (part 3)
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     stress stensor added t.forester may 1994
c     charge group version - t.forester june 1996
c     
c     wl
c     2000/01/18 14:05:36
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      dimension xdf(mxxdf),ydf(mxxdf),zdf(mxxdf)
      dimension fxx(mxatms),fyy(mxatms),fzz(mxatms)
      dimension nexatm(msatms),lexatm(msatms,mxexcl)
      dimension chge(mxatms)
      dimension stress(9)
      
      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/
      data rr3/0.333333333333d0/,r10/0.1d0/,r42/0.02380952381d0/
      data r216/4.62962962963d-3/

#ifdef VAMPIR
      call VTBEGIN(104, ierr)
#endif
c     
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     
c     start of primary loop for forces evaluation
      
      k = 0
      m = -1
      do mm=1,nexatm(ilst),2

        m = m+2
c     
c     atomic index and charge product
        
        iatm=lexatm(ilst,m)
        jatm=lexatm(ilst,m+1)
        chgprd=chge(iatm)*chge(jatm)/epsq*r4pie0
c     
c     calculate interatomic distance

        k = k+1
        rsq=xdf(k)**2+ydf(k)**2+zdf(k)**2
        
        rrr = sqrt(rsq)
        alpr=rrr*alpha
        alpr2=alpr*alpr

c
c     calculate error function and derivative

        if(alpr.lt.1.d-2)then

          erfr=2.d0*chgprd*(alpha/sqrpi)*
     x         (1.d0+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

          egamma=-4.d0*chgprd*(alpha**3/sqrpi)*
     x         (rr3+alpr2*(-2.d0*r10+alpr2*(3.d0*r42-4.d0*alpr2*r216)))

        else

          tt=1.d0/(1.d0+pp*alpha*rrr)
          exp1=exp(-(alpha*rrr)**2)
          erfr =(1.d0-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)*
     x         chgprd/rrr
          egamma = -(erfr-2.d0*chgprd*(alpha/sqrpi)*exp1)/rsq
          
        endif
c     
c     calculate potential energy and virial
        
        engcpe = engcpe-erfr
        vircpe = vircpe - egamma*rsq
c     
c     calculate forces
        
        fx = egamma*xdf(k)
        fy = egamma*ydf(k)
        fz = egamma*zdf(k)
        
        fxx(iatm) = fxx(iatm) + fx
        fyy(iatm) = fyy(iatm) + fy
        fzz(iatm) = fzz(iatm) + fz
        
        fxx(jatm) = fxx(jatm) - fx
        fyy(jatm) = fyy(jatm) - fy
        fzz(jatm) = fzz(jatm) - fz
#ifdef STRESS
c     
c     calculate stress tensor
        
        stress(1) = stress(1)+ xdf(k)*fx
        stress(2) = stress(2)+ xdf(k)*fy
        stress(3) = stress(3)+ xdf(k)*fz
        
        stress(5) = stress(5)+ ydf(k)*fy
        stress(6) = stress(6)+ ydf(k)*fz
        
        stress(9) = stress(9)+ zdf(k)*fz
#endif        
      enddo
#ifdef STRESS
c     
c     complete stress tensor
      
      stress(4) = stress(2)
      stress(7) = stress(3)
      stress(8) = stress(6)
#endif              
      
#ifdef VAMPIR
      call VTEND(104, ierr)
#endif
      return
      end
