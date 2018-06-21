      function thole(a,alp1,alp2,r,n)
      implicit real*8 (a-h,o-z)
c
c     thole's attenuation function
c     exponential decay

      b = (alp1*alp2)**(1.d0/6.d0)

      d = a*r/b

      s = 1.d0
      f = 1.d0
      do i=1,n

        f = f*dble(i)
        s = s + d**i/f

      enddo

      thole = s*dexp(-d)

      return
      end
