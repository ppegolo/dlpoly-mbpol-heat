      subroutine strip(string,length)

c***********************************************************************
c     
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c     
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c     
c     wl
c     1996/02/15 14:33:22
c     1.1.1.1
c     Exp
c     
c***********************************************************************

      character*(*) string
      
      imax = min(length,255)
      do i = 1,imax

        if(string(1:1).eq.' ') then

          do j = 1,imax-1

            string(j:j) = string(j+1:j+1)

          enddo

          string(imax:imax) = ' '

        endif

      enddo

      return
      end
