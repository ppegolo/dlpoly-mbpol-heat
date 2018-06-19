      subroutine forgen
     x  (ltable,idnode,ntpvdw,dlrpot,rcut,ltpvdw,prmvdw,vvv,ggg)
c     
c***********************************************************************
c     
c     dl_poly subroutine for generating potential energy and 
c     force arrays for van der waals forces only
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     modified for t3d dec 1994
c     
c     modifed to be tabulated in r2 space not r space
c     author    - t. forester  march 1993
c     
c     wl
c     2001/05/30 12:40:06
c     1.2
c     Exp
c     
c***********************************************************************
c     
      
#include "dl_params.inc"
      
      logical ltable
      
      dimension ltpvdw(mxvdw),prmvdw(mxvdw,mxpvdw)
      dimension vvv(mxgrid,mxvdw),ggg(mxgrid,mxvdw)
      
c     
c     12 - 6 potential
      
      vv1(r,a,b)=(a/r**6-b)/r**6      
      gg1(r,a,b)=6.d0*(2.d0*a/r**6-b)/r**6
c     
c     lennard-jones potential
      
      vv2(r,a,b)=4.d0*a*(b/r)**6*((b/r)**6-1.d0)
      gg2(r,a,b)=24.d0*a*(b/r)**6*(2.d0*(b/r)**6-1.d0)
c     
c     n - m potential
      
      vv3(r,a,b,c,d)=a/(b-c)*(c*(d/r)**b-b*(d/r)**c)
      gg3(r,a,b,c,d)=a*c*b/(b-c)*((d/r)**b-(d/r)**c)
c     
c     buckingham exp - 6 potential
      
      vv4(r,a,b,c)=a*exp(-r/b)-c/r**6
      gg4(r,a,b,c)=r*a*exp(-r/b)/b-6.d0*c/r**6
c     
c     born-huggins-meyer exp - 6 - 8 potential
      
      vv5(r,a,b,c,d,e)=a*exp(b*(c-r))-d/r**6-e/r**8
      gg5(r,a,b,c,d,e)=r*a*b*exp(b*(c-r))-6.d0*d/r**6-8.d0*e/r**8
      
c     
c     Hydrogen-bond 12 - 10 potential
      
      vv6(r,a,b) = a/r**12 - b/r**10
      gg6(r,a,b) = 12.0d0*a/r**12 - 10.d0*b/r**10
c
c     shifted and force corrected n - m potential (w. smith)
      
      vv7(r,a,b,c,d,b1,c1)=a/(b-c)*(c*(b1**b)*((d/r)**b-
     x  (1.d0/c1)**b)- b*(b1**c)*((d/r)**c-(1.d0/c1)**c)+
     x  b*c*((r/(c1*d)-1.d0)*((b1/c1)**b-(b1/c1)**c)))
      gg7(r,a,b,c,d,b1,c1)=a*c*b/(b-c)*((b1**b)*(d/r)**b-
     x  (b1**c)*(d/r)**c - r/(c1*d)*((b1/c1)**b-(b1/c1)**c))

         
c     morse potential
      vv8(r,a,b,c)=a*((1.d0-exp(-c*(r-b)))**2-1.d0)
      gg8(r,a,b,c)=-2.d0*r*a*c*(1.d0-exp(-c*(r-b)))*exp(-c*(r-b))
c     
c     ttm2 water 12-10-6 potential burnham, jcp116(2002)5115
c     add short range exponential potential JCC
      
!     vv9(r,a,b,c)=a/r**12+b/r**10+c/r**6
!     gg9(r,a,b,c)=12.0d0*a/r**12+10.0d0*b/r**10+6.0d0*c/r**6
      vv9(r,a,b,c,d,e)=a/r**12+b/r**10+c/r**6+d*exp(-e*r)
      gg9(r,a,b,c,d,e)=12.0d0*a/r**12+10.0d0*b/r**10+6.0d0*c/r**6
     $     +d*e*r*exp(-e*r)
!
!     ljex type:
!     Aexp(-Br)+C/r**4+D/r**6
      vv10(r,a,b,c,d)=a*exp(-b*r)+c/r**4+d/r**6
      gg10(r,a,b,c,d)=a*b*r*exp(-b*r)+4.0d0*c/r**4+6.0d0*d/r**6
      
      
#ifdef VAMPIR
      call VTBEGIN(141, ierr)
#endif
c     
c     define grid resolution for potential arrays
      
      dlrpot=rcut**2/dble(mxgrid-4)
      
c     
c     construct arrays for all types of short ranged  potential
      
      do ivdw=1,ntpvdw
        
        if(ltpvdw(ivdw).eq.1)then
          
          do i=1,mxgrid
            
            rrr=sqrt(dble(i)*dlrpot)
            vvv(i,ivdw)=vv1(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
            ggg(i,ivdw)=gg1(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
            
          enddo
          
        else if(ltpvdw(ivdw).eq.2)then
          
          do i=1,mxgrid
            
            rrr=sqrt(dble(i)*dlrpot)
            vvv(i,ivdw)=vv2(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
            ggg(i,ivdw)=gg2(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
            
          enddo
          
        else if(ltpvdw(ivdw).eq.3)then

          an = prmvdw(ivdw,2)
          am = prmvdw(ivdw,3)
          if(an.le.am) call error(idnode,470)
          
          do i=1,mxgrid
            
            rrr=sqrt(dble(i)*dlrpot)
            vvv(i,ivdw)=vv3(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x        prmvdw(ivdw,3),prmvdw(ivdw,4))
            ggg(i,ivdw)=gg3(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x        prmvdw(ivdw,3),prmvdw(ivdw,4))
            
          enddo
          
        else if(ltpvdw(ivdw).eq.4)then
          
          do i=1,mxgrid
            
            rrr=sqrt(dble(i)*dlrpot)
            vvv(i,ivdw)=vv4(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x        prmvdw(ivdw,3))
            ggg(i,ivdw)=gg4(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x        prmvdw(ivdw,3))
            
          enddo
          
        else if(ltpvdw(ivdw).eq.5)then
          
          do i=1,mxgrid
            
            rrr=sqrt(dble(i)*dlrpot)
            vvv(i,ivdw)=vv5(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x        prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5))
            ggg(i,ivdw)=gg5(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x        prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5))
            
          enddo
          
        else if(ltpvdw(ivdw).eq.6) then
          
          do i = 1,mxgrid
            
            rrr=sqrt(dble(i)*dlrpot)
            vvv(i,ivdw)=vv6(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
            ggg(i,ivdw)=gg6(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2))
            
          enddo

        else if(ltpvdw(ivdw).eq.7) then
          
          r0 = prmvdw(ivdw,4)
          an = prmvdw(ivdw,2)
          am = prmvdw(ivdw,3)
          eps = prmvdw(ivdw,1)
 
          if(an.le.am) call error(idnode,470)

          gamma = rcut/r0
          if(gamma.lt.1.d0) call error(idnode,468)
          beta= gamma*((gamma**(am+1.d0)-1.d0)/(gamma**(an+1.d0)-1.d0))
     x      **(1.d0/(an-am))
          alpha = -(an-am)/(am*(beta**an)*(1.d0+(an/gamma-an-1.d0)
     x      /gamma**an) -   an*(beta**am)*(1.d0+(am/gamma-am-1.d0)
     x      /gamma**am))
          eps = eps*alpha

          do i = 1,mxgrid
            
            rrr=sqrt(dble(i)*dlrpot)
            vvv(i,ivdw)=vv7(rrr,eps,an,am,r0,beta,gamma)
            ggg(i,ivdw)=gg7(rrr,eps,an,am,r0,beta,gamma)

          enddo
          
       else if(ltpvdw(ivdw).eq.8) then
          
          do i = 1,mxgrid
             
             rrr=dble(i)*dlrpot
             vvv(i,ivdw)=vv8(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3))
             ggg(i,ivdw)=gg8(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3))
             
          enddo
          
       else if(ltpvdw(ivdw).eq.9) then
          
          do i=1,mxgrid
             
c$$$  prmvdw(ivdw,1) = -1.329565985d+6 * 4.184d2
c$$$  prmvdw(ivdw,2) =  3.632560798d+5 * 4.184d2
c$$$  prmvdw(ivdw,3) = -2.147141323d+3 * 4.184d2
             
             rrr=dble(i)*dlrpot
             vvv(i,ivdw)=vv9(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5))
             ggg(i,ivdw)=gg9(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5))
             
          enddo
          
       else if(ltpvdw(ivdw).eq.10) then
          
          do i=1,mxgrid
             
             rrr=dble(i)*dlrpot
             vvv(i,ivdw)=vv10(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3),prmvdw(ivdw,4))
             ggg(i,ivdw)=gg10(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3),prmvdw(ivdw,4))
             
          enddo
          
       else if(ltpvdw(ivdw).eq.11) then
          
          do i=1,mxgrid
             
             rrr=dble(i)*dlrpot
             vvv(i,ivdw)=vv11(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5),prmvdw(ivdw,6))
             ggg(i,ivdw)=gg11(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3),prmvdw(ivdw,4),prmvdw(ivdw,5),prmvdw(ivdw,6))
             
          enddo
          
       else if(ltpvdw(ivdw).eq.12)then
          
          do i=1,mxgrid
             
             rrr=dble(i)*dlrpot
             vvv(i,ivdw)=vv12(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3))
             ggg(i,ivdw)=gg12(rrr,prmvdw(ivdw,1),prmvdw(ivdw,2),
     x            prmvdw(ivdw,3))
             
          enddo
          
       else if(ltpvdw(ivdw).lt.100) then
          
          if(.not.ltable)call error(idnode,150)
          
        endif
        
      enddo
      
#ifdef VAMPIR
      call VTEND(141, ierr)
#endif
      return

      contains
!
!     TTM4-F dispersion potential
!

      real(8) function vv11(r,A6,A8,A10,A12,A14,A16)

      implicit none

      real(8) :: r,A6,A8,A10,A12,A14,A16

      real(8) :: x

      x = 1.d0/r/r

      vv11 =
     x    x*x*x*(A6 + x*(A8 + x*(A10 + x*(A12 + x*(A14 + x*A16)))))

      end function vv11

      ! -r*(d/dr)
      ! v11(r)

      real(8) function gg11(r,A6,A8,A10,A12,A14,A16)

      implicit none

      real(8) :: r,A6,A8,A10,A12,A14,A16
      real(8) :: x

      x = 1.d0/r/r

      gg11 = x*x*x*(6.d0*A6 + x*(8.d0*A8 + x*(10.d0*A10
     x      + x*(12.d0*A12 + x*(14.d0*A14 + x*16.d0*A16)))))

      end function gg11

!
!     Modified Buckingham for two-body: ttmx
! 

      real(8) function vv12(r,a,b,c)
      implicit none
      real(8) :: r,a,b,c

      vv12=a/r*exp(-r/b)-c/r**6

      end function vv12

      ! -r*(d/dr) v12(r)

      real(8) function gg12(r,a,b,c)
      implicit none
      real(8) :: r,a,b,c

      gg12=a*exp(-r/b)*(1.d0/b+1.d0/r)-6.d0*c/r**6

      end function gg12

      end
