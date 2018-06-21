      subroutine scl_csum(nnn,tot,aaa)

c**********************************************************************
c
c     dl_poly subroutine to calculate the scalar sum of the elements
c     of a complex array
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      complex*16 aaa(*),tot

      tot=(0.d0,0.d0)

      do i=1,nnn

        tot=tot+aaa(i)

      enddo

      return
      end
