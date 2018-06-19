      subroutine cpy_rtc(nnn,aaa,bbb)

c**********************************************************************
c
c     dl_poly subroutine for copying a real array into a complex array
c     of the same dimension
c
c     copyright daresbury laboratory 1998
c     author w.smith oct 1998
c
c**********************************************************************

      real*8 aaa(*)
      complex*16 bbb(*)

      do i=1,nnn

        bbb(i)=cmplx(aaa(i),0.d0,kind=8)

      enddo

      return
      end

