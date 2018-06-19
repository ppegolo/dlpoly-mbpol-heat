      subroutine gstate(check)
c     
c***********************************************************************
c     
c     dl_poly global status subroutine : gisum version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     MPI version -  t. forester may 1995
c     
c     wl
c     1.2
c     Exp
c***********************************************************************
c     

      logical check

#ifndef SERIAL

      i = 0
      if(.not.check) i = 1

      call gisum(i,1,idum)
      
      check = (i.eq.0)

#endif
      
      return
      end
