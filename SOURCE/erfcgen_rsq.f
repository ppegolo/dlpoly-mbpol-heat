      subroutine erfcgen
     x  (keyfce,alpha,drewd,rcut,erc,fer)

c***********************************************************************
c     
c     dlpoly routine for generating interpolation tables for 
c     erfc and its derivative - for use with ewald sum.
c     r-squared tabulation
c     
c     copyright daresbury laboratory 1994
c     author t.forester dec 1994
c     
c     wl
c     2001/06/12 12:44:43
c     1.4
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      dimension erc(mxegrd),fer(mxegrd)
      
      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/
#ifdef VAMPIR
      call VTBEGIN(136, ierr)
#endif
c     
c     look-up tables for real space part of ewald sum
      
      if(keyfce/2.eq.1.or.keyfce/2.eq.6) then
        
        drewd = rcut**2/dble(mxegrd-4)
        
        do i = 1,mxegrd
          
          rrr = sqrt(dble(i)*drewd)
          rsq = rrr*rrr
          tt=1.d0/(1.d0+pp*alpha*rrr)
          exp1=exp(-(alpha*rrr)**2)
          erc(i)=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))
     x      *exp1/rrr
          fer(i)=erc(i)/rsq + 2.d0*(alpha/sqrpi)*exp1/rsq
          
        enddo
        
      endif

#ifdef VAMPIR
      call VTEND(136, ierr)
#endif
      return
      end
