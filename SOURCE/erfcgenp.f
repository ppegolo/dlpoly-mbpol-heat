      subroutine erfcgenp
     x  (keyfce,alpha,drewd,rcut,ercp)

c***********************************************************************
c     
c     dlpoly routine for generating interpolation tables for 
c     erfc and its derivative - for use with ewald sum.
c     
c     copyright daresbury laboratory 1994
c     author t.forester dec 1994
c     
c     wl
c     2001/06/12 12:43:54
c     1.4
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      dimension ercp(mxegrd,0:3)
      
#ifdef VAMPIR
      call VTBEGIN(135, ierr)
#endif
c     
c     look-up tables for real space part of ewald sum
      
      if(keyfce/2.eq.1.or.keyfce/2.eq.6) then
        
        drewd = rcut/dble(mxegrd-4)
        
        do i = 1,mxegrd
          
          rrr = dble(i)*drewd
          rsqi = 1.d0/rrr/rrr 
          a=alpha*rrr
          ercp(i,0)=(erfcc(a)-1.d0)/rrr
          exp2a=exp(-a*a)

          do n=1,3

             fn=dble(n)
             ercp(i,n)=rsqi*((2.d0*fn-1.d0)*ercp(i,n-1)
     x          +((2.d0*alpha**2)**fn)*exp2a/alpha/sqrpi)

          enddo

        enddo
        
      endif

#ifdef VAMPIR
      call VTEND(135, ierr)
#endif
      return
      end
