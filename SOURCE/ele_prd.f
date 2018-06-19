      subroutine ele_prd(nnn,aaa,bbb,ccc)

c**********************************************************************
c
c     dl_poly subroutine for element by element product of
c     a real array (bbb) and a complex array (ccc)
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      real*8 bbb(*)
      complex*16 aaa(*),ccc(*)

      do i=1,nnn

        aaa(i)=bbb(i)*ccc(i)

      enddo

      return
      end
