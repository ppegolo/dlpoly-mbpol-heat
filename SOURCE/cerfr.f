      function cerfr(alpha,rrr)

c***********************************************************************
c     
c     dl_poly function for generating complementary error function
c     divided by r
c     
c     copyright - daresbury laboratory 2001
c     author    - w. smith february 2001
c     
c     wl
c     2001/06/12 13:10:08
c     1.1
c     Exp
c     
c***********************************************************************

      implicit real*8(a-h,o-z)

      sqrpi=1.7724538509055159d0

c     starting values

      h0=0.d0
      h1=(2.d0*alpha/sqrpi)*exp(-(alpha*rrr)**2)

c     generate function by recursion

      rr2=rrr*rrr
      do k=300,1,-1

        fac=(2.d0*alpha**2)**(k-1)
        h0=(h0*rr2+fac*h1)/dble(2*k-1)

      enddo
        
      cerfr=1.d0/rrr-h0

      return
      end
