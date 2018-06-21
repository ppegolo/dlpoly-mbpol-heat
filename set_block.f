      subroutine set_block(nnn,ccc,aaa)

c**********************************************************************
c
c     dl_poly subroutine to initialise an array to a single value
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      real*8 ccc,aaa(nnn)

      do i=1,nnn,2

        aaa(i)=ccc
        aaa(i+1)=ccc

      enddo
      
      return
      end
